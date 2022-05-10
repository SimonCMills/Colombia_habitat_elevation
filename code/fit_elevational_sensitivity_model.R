# fitting script 

library(dplyr); library(flocker); library(brms); library(ape)

# read data
c_phylo <- readRDS("output/consensus_phylogeny.rds")
adf <- readRDS("output/analysis_dataset.rds")

sp_list <- unique(adf$species_eltontraits)

# generate distance matrix
drop_index <- which(!(c_phylo$tip.label %in% adf$species_eltontraits))
phylo_subset <- drop.tip(c_phylo, drop_index)
ultratree <- chronopl(phylo_subset, lambda = 0, age.min = 1)

A_1 <- vcv(phylo_subset)
A <- A_1/max(A_1)

det <- adf %>% select(d1:d4) %>% 
    mutate_all(function(x)ifelse(x==-99, NA, x)) %>%
    as.matrix

v_cov <- list(time = adf %>% select(v1:v4) %>% 
                  mutate_all(function(x)ifelse(x==-99, NA, x)) %>%
                  as.matrix, 
              observer = adf %>% select(o1: o4) %>% 
                  mutate_all(function(x)ifelse(x==-99, NA, x)) %>%
                  as.matrix) 

cu_cov <- adf %>%
    mutate(range_pos_upr = ifelse(range_pos > 0, range_pos, 0), 
           medium_dependency = ifelse(forest_dep == "Medium", -1, 1),
           range_size_sc = scale(scale_both)) %>%
    select(rotation, 
           species, 
           species_eltontraits,
           id_cl_sp, 
           id_site_sp, 
           habitat_std, 
           range_pos, 
           range_pos2,
           range_size_sc,
           range_pos_upr,
           ele_pt_std,
           ele_midpoint_std, 
           obsvr_point, 
           id_obs_sp, 
           medium_dependency)

fd1 <- make_flocker_data(det, cu_cov, v_cov)

fit_prior <- c(brms::set_prior("normal(-2, 1.5)", class = "Intercept", dpar="occ"), 
               brms::set_prior("normal(0, .5)", class = "b", dpar = "occ"), # dependency interactions
               brms::set_prior("normal(0, 2)", coef = "range_pos", dpar = "occ"), 
               brms::set_prior("normal(-4, 2)", coef = "range_pos2", dpar = "occ"), 
               brms::set_prior("normal(0, 2)", coef = "ele_midpoint_std", dpar = "occ"), 
               brms::set_prior("normal(0, 2)", coef = "ele_midpoint_std:habitat_std", dpar = "occ"), 
               brms::set_prior("normal(0, 2)", coef = "range_pos:habitat_std", dpar = "occ"), 
               brms::set_prior("normal(0, 2)", coef = "habitat_std:range_pos_upr", dpar = "occ"), 
               brms::set_prior("normal(0, 1)", class = "b", dpar=""), 
               brms::set_prior("normal(-2, 2)", class = "Intercept", dpar = ""), 
               brms::set_prior("normal(0, 1)", class = "sd"),
               brms::set_prior("normal(0, 3)", class = "sd", group = "id_cl_sp", dpar = "occ"),
               brms::set_prior("normal(0, 3)", class = "sd", group = "id_site_sp", dpar = "occ"))

## fit model 
dirname <- paste0(getwd(), "/stan_out/elev_sens_mod1")
if(!dir.exists(dirname)) dir.create(dirname)

writename <- paste0(format(Sys.time(), "warmup_%d-%m-%Y_%H%M_j"), job_id)

fit6 <- flocker::flock(f_occ = ~ 1 + range_pos + range_pos2 + ele_midpoint_std + 
                           habitat_std + ele_midpoint_std:habitat_std + 
                           range_pos:habitat_std +
                           range_pos_upr:habitat_std +
                           # interaction terms for medium_dependency species
                           medium_dependency + medium_dependency:range_pos +
                           medium_dependency:range_pos2 + 
                           medium_dependency:ele_midpoint_std + 
                           medium_dependency:habitat_std + 
                           medium_dependency:habitat_std:ele_midpoint_std + 
                           medium_dependency:range_pos:habitat_std +
                           medium_dependency:range_pos_upr:habitat_std +
                           range_size_sc + range_size_sc:habitat_std +
                           # random effects
                           (1 + habitat_std + range_pos + range_pos2||species) + 
                           (1|id_cl_sp) + (1|id_site_sp) + 
                           # phylogenetic terms
                           (1 + habitat_std|gr(species_eltontraits, cov=A)), 
                       f_det = ~ 1 + time + (1 + time|species) + obsvr_point + 
                           (1|rotation) + (1|id_obs_sp), 
                       data2 = list(A = A),
                       flocker_data = fd1, 
                       rep_constant = F, 
                       prior = fit_prior, 
                       # run info
                       backend = "cmdstanr", 
                       chains = 4, cores = 4,
                       iter = 2000, 
                       max_treedepth = 9, 
                       adapt_delta = .95,
                       # save info
                       file = paste0(dirname, "/", writename), 
                       save_warmup = T,
                       output_dir = dirname, 
                       output_basename = writename)
