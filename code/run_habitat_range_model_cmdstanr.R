## run elevation model

## Housekeeping ----
library("dplyr"); library("ggplot2"); library("reshape2"); library("bayesplot")
library("cmdstanr"); library("posterior")

# get cores/threads from bash
cpu_info <- as.numeric(commandArgs(T))
n_cores <- n_chains <- cpu_info[1]
n_threads <- cpu_info[2]

####
range_scaling <- function(ele, mid, scale) {
    (ele - mid)/(scale)
}
range_unscaling <- function(ele_scaled, mid, scale) {
    ele_scaled *scale + mid
}


## Read in datasets ----
df_bird <- readRDS("data/analysisDataset_EasternCordillera.rds") %>%
    filter(Visit != 5)

#* forest and landscape metrics ----
forest_cover <- readRDS("data/hansen_500m_buffer.rds") %>%
    mutate(tc_2000 = treecover2000 > 50, 
           tc_2018 = as.numeric(tc_2000 - (lossyear != 0) == 1)) %>%
    group_by(point_id) %>%
    summarise(amt_tc = sum(tc_2018), pct_tc = amt_tc/n())

lscape_metrics <- readRDS("data_landscapeMetrics_draft_v1.rds") %>%
    group_by(Point) %>% slice(1)

#* elevational ranges (Quinones) ----
ele_range <- read.csv("data/elevational_ranges_Quinones.csv") %>%
    mutate(Species = paste0(Genus, "_", Species)) %>%
    mutate(Lower = ifelse(is.na(Lower), 0, Lower))

#* point info ----
df_ptInfo <- read.csv("data/CO_sampling_points_metafile_ver2.csv") %>%
    as_tibble() %>%
    left_join(., forest_cover) %>%
    mutate(hab = substr(point_id, 3, 3)) %>%
    rename(Point = point_id) %>%
    left_join(., lscape_metrics) %>%
    # filter(Point %in% unique_pt) %>%
    dplyr::select(Point, lat, long, ele=ALOSelev, amt_tc, pct_tc, Habitat, 
                  pct_nat_hab = Prop_natural_habitat)

## Get full df and zf version ----
# filtering for forest points
df_det <- df_bird %>% 
    left_join(., df_ptInfo) %>%
    filter(Habitat == "Forest") %>%
    group_by(Point, Visit) %>%
    filter(!duplicated(Species)) %>%
    ungroup %>%
    mutate(Q = 1) %>%
    dplyr::select(Site, Species, Point, Visit, Time, Q) 

# note: 72 point:visit combinations that have duplicated species-- FIX 
# det %>% group_by(Point, Visit) %>%
#     filter(duplicated(Species))

uniqueVisits <- with(df_det, 
                     expand.grid(Point = unique(Point), Visit = 1:4)) %>%
    left_join(., df_bird[c("Point", "Visit", "Time", "Site")] %>% unique) %>%
    mutate(Time = ifelse(is.na(Time), median(Time, na.rm=T), Time))

uniqueSpecies <- unique(df_det$Species)
n_species <- length(uniqueSpecies)

# zero-filled dataframe
zf_det <- replicate(n_species, uniqueVisits, simplify=F) %>%
    bind_rows(., .id="Species") %>%
    mutate(Species = uniqueSpecies[as.numeric(Species)]) %>%
    anti_join(., df_det) %>%
    mutate(Q = 0) %>%
    bind_rows(., df_det)

#* Take subset ----
# now filter with species on at least 5 points
nPt_bySpecies <- zf_det %>% 
    filter(Q == 1) %>%
    dplyr::select(Species, Point) %>%
    group_by(Species) %>%
    summarise(n_pt = length(unique(Point))) %>% 
    arrange(desc(n_pt)) %>%
    filter(n_pt >= 20)

zf_subset <- filter(zf_det, Species %in% nPt_bySpecies$Species) %>%
    mutate(Time_sc = (Time-median(Time))/sd(Time))

det <- dcast(zf_subset, Point + Species ~ Visit, value.var="Q")
vis <- dcast(zf_subset, Point + Species ~ Visit, value.var="Time_sc")

eles <- left_join(det, df_ptInfo) %>%
    left_join(., ele_range) %>%
    mutate(midpoint = (Lower + Upper)/2, 
           scale = (Upper - Lower)/2, 
           ele_sc = range_scaling(ele, midpoint, scale), 
           id_site_sp = interaction(substr(Point, 1, 2), Species), 
           id_site_sp_num = as.numeric(id_site_sp))

# Deal with cases where there is missing Upper Lower info
eles2 <- eles %>% filter(!is.na(eles$ele_sc))
det2 <- det %>% filter(!is.na(eles$ele_sc))
vis2 <- vis %>% filter(!is.na(eles$ele_sc))

Q <- as.numeric(rowSums(det2[,-c(1,2)]) >= 1)

ele_forPred <- eles2 %>% group_by(Species) %>% 
    summarise(min_ele = min(ele_sc), max_ele=max(ele_sc)) %>%
    apply(., 1, function(x) data.frame(ele_sc = seq(x["min_ele"], x["max_ele"], len=50))) %>%
    bind_cols

lookup <- tibble(Species = as.factor(det2$Species), 
                 id_sp = as.numeric(as.factor(det2$Species))) %>% unique %>%
    mutate(id_sp = as.numeric(id_sp))

stan_data <- list(n_visit = 4,
                  n_species = length(unique(det2$Species)),
                  n_points = length(unique(det2$Point)),
                  n_preds = 50,
                  n_tot = nrow(det2),
                  id_sp = as.numeric(as.factor(det2$Species)),
                  id_pt = as.numeric(as.factor(det2$Point)),
                  id_site_sp = eles2$id_site_sp_num,
                  n_sites_species = max(eles2$id_site_sp_num),
                  Q = Q,
                  vis_cov1 = vis2[,-c(1,2)],
                  det_data = det2[,-c(1,2)],
                  ele = eles2$ele_sc,
                  ele_forPred = ele_forPred,
                  habitat = 1 - eles2$pct_tc,
                  grainsize=1)

## Run mod ----
mod <- cmdstan_model("code/stan_files/elevation_detection_habitat.stan", 
                     cpp_options = list(stan_threads = T))

samps <- mod$sample(data = stan_data, 
                    chains = 3, 
                    parallel_chains = n_chains, #n_chains, 
                    threads_per_chain = n_threads, #n_threads, 
                    iter_warmup = 1000, 
                    iter_sampling = 2000)

draws <- samps$draws()
# draws without junk
draws_clean_1 <- draws[,,!grepl("raw|ele_scaled", variables(draws))]

# draws without junk or site effects -- retain
saveRDS(samps, "output/model_habitatModel.rds")
saveRDS(draws_clean_1, "output/draws_habitatModel.rds")
