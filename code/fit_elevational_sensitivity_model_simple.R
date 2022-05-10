# script to fit model to raw data and plot for SOM

library(brms); library(dplyr); library(ggplot2)

fname <- "output/model_without_detection_correction.rds"
run_mod <- !file.exists(fname)

if(run_mod) {
    fit_prior <- c(brms::set_prior("normal(-2.5, 1.5)", class = "Intercept"), 
                   brms::set_prior("normal(0, .5)", class = "b"), # dependency interactions
                   brms::set_prior("normal(0, 2)", coef = "range_pos"), 
                   brms::set_prior("normal(-4, 2)", coef = "range_pos2"), 
                   brms::set_prior("normal(0, 2)", coef = "ele_midpoint_std"), 
                   brms::set_prior("normal(0, 2)", coef = "ele_midpoint_std:habitat_std"), 
                   brms::set_prior("normal(0, 2)", coef = "habitat_std:range_pos"), 
                   brms::set_prior("normal(0, 2)", coef = "habitat_std:range_pos_upr"), 
                   brms::set_prior("normal(0, 1)", class = "sd"),
                   brms::set_prior("normal(0, 3)", class = "sd", group = "id_cl_sp"),
                   brms::set_prior("normal(0, 3)", class = "sd", group = "id_site_sp"), 
                   brms::set_prior("normal(0, 1)", class = "sd", group = "species_eltontraits"))
    
    f2 <- brm(Q ~ ele_midpoint_std + habitat_std + 
                  ele_midpoint_std:habitat_std + 
                  range_pos + range_pos2 + 
                  range_pos:habitat_std + 
                  range_pos_upr:habitat_std + 
                  medium_dependency + 
                  medium_dependency:habitat_std +
                  medium_dependency:habitat_std:ele_midpoint_std + 
                  medium_dependency:range_pos +
                  medium_dependency:range_pos2 +
                  medium_dependency:range_pos:habitat_std +
                  medium_dependency:range_pos_upr:habitat_std + 
                  # random effects
                  (1 + habitat_std + range_pos + range_pos2|species) + 
                  (1|id_cl_sp) + (1|id_site_sp) + 
                  # phylogenetic terms
                  (1 + habitat_std|gr(species_eltontraits, cov=A)), 
              prior = fit_prior,
              adapt_delta = .97, 
              data = adf2, family="bernoulli",
              file = "output/model_without_detection_correction",
              data2 = list(A = A), 
              backend="cmdstanr", chains=4, cores=4)
} else {
    f2 <- readRDS("output/model_without_detection_correction.rds")
}

# extract coefficients & plot
# note: this is just tweaked version of code from plotting script
feffs <- fixef(f2, summary=F) %>%
    as_tibble %>%
    slice(seq(1, 4000, 8))

names(feffs)
high_dep_effs <- feffs %>%
    as_tibble() %>%
    mutate(inter_eff = `ele_midpoint_std:habitat_std`, 
           intra_eff = `habitat_std:range_pos`, 
           intra_eff_upper = `habitat_std:range_pos` + `habitat_std:range_pos_upr`) %>%
    summarise_at(vars(inter_eff, intra_eff, intra_eff_upper), 
                 list(mean = mean, 
                      lwr=function(x) quantile(x, .025), 
                      upr=function(x) quantile(x, .975), 
                      pd = function(x) round(min(sum(x<0), sum(x>0))/length(x), 2)
                 )) %>%
    mutate(dep = "high") %>%
    reshape2::melt(., id.vars="dep") %>%
    mutate(type = gsub(".*_(.*$)", "\\1", variable), 
           par = gsub("(.*)_.*$", "\\1", variable)) %>%
    reshape2::dcast(., formula = dep + par ~ type)

med_dep_effs <- feffs %>%
    as_tibble() %>%
    mutate(inter_eff = `ele_midpoint_std:habitat_std` + `ele_midpoint_std:habitat_std:medium_dependency`, 
           intra_eff = `habitat_std:range_pos` + `habitat_std:range_pos:medium_dependency`, 
           intra_eff_upper = intra_eff + `habitat_std:range_pos_upr` + 
               `habitat_std:medium_dependency:range_pos_upr`) %>%
    summarise_at(vars(inter_eff, intra_eff, intra_eff_upper), 
                 list(mean = mean, 
                      lwr=function(x) quantile(x, .025), 
                      upr=function(x) quantile(x, .975), 
                      pd = function(x) round(min(sum(x<0), sum(x>0))/length(x),2)
                 )) %>%
    mutate(dep = "medium") %>%
    reshape2::melt(., id.vars="dep") %>%
    mutate(type = gsub(".*_(.*$)", "\\1", variable), 
           par = gsub("(.*)_.*$", "\\1", variable)) %>%
    reshape2::dcast(., formula = dep + par ~ type)


ylabs <- expression(
    beta[8]~" (elevational midpoint x habitat)",
    beta[7]~" (scaled elevation[lower] x habitat)",
    beta[6]~" (scaled elevation[upper] x habitat)",
)

effs_p1 <- bind_rows(high_dep_effs, med_dep_effs) %>%
    mutate(pd = ifelse(pd  == 0, "<0.01", pd)) %>%
    ggplot(aes(y = mean, ymin=lwr, ymax=upr, x=par, fill = dep, group=dep)) + 
    geom_errorbar(width=0, position = position_dodge(width=.5)) +
    geom_point(position = position_dodge(width=.5), pch=21, size=2) +
    coord_flip()  +
    theme_bw() +
    scale_fill_manual(values=c("black", "white")) +
    scale_x_discrete(labels = ylabs) +
    scale_y_continuous(expand=c(.1,.1)) +
    theme(axis.text = element_text(colour="black"), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"),
          # panel.background = element_rect(fill="grey90"),
          aspect.ratio = .7,
          panel.grid = element_blank(), 
          legend.position=c(.9,.89), 
          legend.background = element_rect(colour="black"), 
          legend.key.height=unit(.9,"line"), 
          plot.margin = unit(c(0.2, 2, 0.2, 0.2), units = "cm")
    ) +
    #guides(fill="none") +
    geom_hline(yintercept = 0, lty="longdash") +
    labs(y = "", x = "", fill = "Dependency") +
    geom_text(aes(label=pd, y = upr), position = position_dodge(width=.5), hjust=-0.2, vjust=.3, size=3)


# summ_fit <- summary(fit)
sigmas <- VarCorr(f2, summary=F)
lambdas <- (sigmas$species_eltontraits$sd/
                (sigmas$species$sd[,3:4] + sigmas$species_eltontraits$sd))[seq(1, 4000, 8),] %>%
    as_tibble

# lambdas
lambdas2 <- lambdas %>% 
    summarise_all(list(mean = mean, lwr = function(x)quantile(x, .025), 
                       upr = function(x) quantile(x, .975))) %>%
    reshape2::melt(.) %>%
    mutate(type = gsub(".*_(.*$)", "\\1", variable), 
           par = gsub("(.*)_.*$", "\\1", variable))  %>%
    reshape2::dcast(., formula = par ~ type)

effs_p2 <- lambdas2 %>%
    filter(par == "habitat_std") %>%
    ggplot(aes(x = par, y = mean,  ymin=lwr, ymax=upr)) +
    geom_errorbar(width=0, position = position_dodge(width=.5)) +
    geom_point(position = position_dodge(width=.5), pch=21, size=2, fill="black") +
    coord_flip()  +
    theme_bw() +
    scale_fill_manual(values=c("black", "white")) +
    scale_x_discrete(labels = expression(lambda[habitat])) +
    scale_y_continuous(expand=c(.01,.01), limits=c(0,1), breaks=seq(0, 1, .1)) +
    theme(axis.text = element_text(colour="black"), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"),
          # panel.background = element_rect(fill="grey90"),
          panel.grid = element_blank(), 
          legend.position=c(.86,.15), 
          legend.background = element_rect(colour="black"), 
          legend.key.height=unit(.9,"line"), 
          plot.margin = unit(c(0.2, 2, 0.2, 0.2), units = "cm")
    ) +
    geom_hline(yintercept=c(0,1), lty="longdash") +
    labs(x="", y="Parameter value (\u00B1 90% CI)")

effs_both <- egg::ggarrange(effs_p1, effs_p2, ncol=1, heights=c(1, .1))

ggsave("figures/supplementary_analysis_coefficient_ests.png", plot = effs_both, 
       units = "mm", width=250, height=190, dpi=400)
