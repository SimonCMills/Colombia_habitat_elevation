# plotting script for main figures

library(flocker); library(brms); library(dplyr); library(ggplot2)

# read in model and dataset ----
fit <- readRDS("output/model_elevational_sensitivity.rds")
analysis_df <- readRDS("output/analysis_df.rds")
dep_lvls <- c(high = 1, medium = -1)

# average effects ----
newdat <- expand.grid(range_pos = seq(-1, 1, len=30), 
                      ele_midpoint = round(seq(800, 2700, len=4), -2), 
                      habitat = c(1, .6), 
                      range_size_sc = mean(scale(analysis_df$scale_both)),
                      medium_dependency = dep_lvls, 
                      time = 1, 
                      obsvr_point = "JS") %>%
    mutate(range_pos2 = range_pos^2, 
           range_pos_upr = ifelse(range_pos < 0, 0, range_pos),
           ele_midpoint_std = (ele_midpoint - mean(analysis_df$midpoint_both))/sd(analysis_df$midpoint_both), 
           habitat_std = (habitat - mean(analysis_df$habitat))/sd(analysis_df$habitat), 
           ele_pt = 700 * range_pos + ele_midpoint_std, 
           ele_pt_std = (ele_pt - mean(analysis_df$ele_pt))/sd(analysis_df$ele_pt))

preds <- fitted_flocker(fit, type = "occupancy", 
                       new_data = newdat, re_formula = NA, summarise = F, 
                       ndraws = 500)
## calculate log odds ---- 
logOdds <- cbind(newdat, boot::logit(preds)) %>%
    group_by(ele_midpoint, range_pos, medium_dependency) %>%
    summarise_at(vars(paste0("iter_", 1:500)), .funs = function(x) x[2]-x[1])

logOdds2 <-logOdds %>%
    ungroup %>%
    dplyr::select(iter_1:iter_500) %>%
    as.matrix(.)

logOdds3 <- logOdds %>%
    ungroup %>%
    mutate(mean = rowMeans(logOdds2), 
           lwr = matrixStats::rowQuantiles(logOdds2, probs = .05),
           upr = matrixStats::rowQuantiles(logOdds2, probs = .95)) %>%
    dplyr::select(-(iter_1:iter_500))

## summarise occupancy ----
median_occ <- cbind(newdat, preds)
median_occ_mat <- median_occ %>%
    ungroup %>%
    dplyr::select(iter_1:iter_500) %>%
    as.matrix(.)

median_occ2 <- median_occ %>%
    ungroup %>%
    mutate(mean = rowMeans(median_occ_mat), 
           lwr = matrixStats::rowQuantiles(median_occ_mat, probs = .05),
           upr = matrixStats::rowQuantiles(median_occ_mat, probs = .95)) %>%
    dplyr::select(-(iter_1:iter_500))

## plot ----
cols <- RColorBrewer::brewer.pal(8, name = "RdBu")[c(1,8)]#viridis::viridis(n=4)

occ_p1 <- median_occ2 %>%
    filter(medium_dependency == dep_lvls["high"], !(ele_midpoint == 800)) %>%
    mutate(label = case_when(ele_midpoint == 1400 ~ "(a) high dependency, 1400 m", 
                             ele_midpoint == 2100 ~ "(b) high dependency, 2100 m", 
                             ele_midpoint == 2700 ~ "(c) high dependency, 2700 m")) %>%
    ggplot(., aes(range_pos, mean, group=habitat_std)) +
    geom_line(aes(col=factor(habitat_std))) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, fill=interaction(habitat_std)), alpha=.1) +
    facet_wrap(~label, ncol=4) +
    guides(fill="none", lty="none", col="none") +
    theme_bw() +
    theme(axis.text = element_text(colour="black"), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"),
          # panel.background = element_rect(fill="grey90"),
          aspect.ratio = .7,
          panel.grid = element_blank(), 
          legend.position=c(.955,.37), 
          legend.background = element_blank(), 
          legend.key.height=unit(.9,"line"), 
          axis.text.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult=.03)) +
    scale_x_continuous(expand = expansion(mult=.03)) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    labs(x = "", 
         y = "P(Occupancy)", colour = "", fill= "")

occ_p2 <- median_occ2 %>%
    filter(medium_dependency == dep_lvls["medium"], !(ele_midpoint == 800)) %>%
    mutate(label = case_when(ele_midpoint == 1400 ~ "(d) medium dependency, 1400 m", 
                             ele_midpoint == 2100 ~ "(e) medium dependency, 2100 m", 
                             ele_midpoint == 2700 ~ "(f) medium dependency, 2700 m")) %>%
    ggplot(., aes(range_pos, mean, group=habitat_std)) +
    geom_line(aes(col=factor(habitat_std)), linetype=2) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, fill=interaction(habitat_std)), alpha=.1) +
    facet_wrap(~label, ncol=4) +
    guides(fill="none", lty="none", col="none") +
    theme_bw() +
    theme(axis.text = element_text(colour="black"), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"),
          # panel.background = element_rect(fill="grey90"),
          aspect.ratio = .7,
          panel.grid = element_blank(), 
          legend.position=c(.955,.37), 
          legend.background = element_blank(), 
          legend.key.height=unit(.9,"line"), 
          axis.text.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult=.03)) +
    scale_x_continuous(expand = expansion(mult=.03)) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    labs(x = "", 
         y = "P(Occupancy)", colour = "", fill= "")

occ_p3 <- logOdds3 %>%
    filter(!(ele_midpoint == 800)) %>%
    mutate(label = case_when(ele_midpoint == 1400 ~ "(h) log-odds, 1400 m", 
                             ele_midpoint == 2100 ~ "(i) log-odds, 2100 m", 
                             ele_midpoint == 2700 ~ "(j) log-odds, 2700 m")) %>%
    ggplot(aes(range_pos, mean, group=factor(medium_dependency), 
                     lty = factor(medium_dependency))) + 
    geom_line() +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha = .1) +
    scale_linetype_manual(values=c(2, 1)) +
    # scale_y_continuous(breaks = seq(-4, 4, 2)) +
    facet_wrap(~label, ncol=4) +
    geom_hline(yintercept = 0, lty="longdash") +
    theme_bw() +
    theme(axis.text = element_text(colour="black"), 
          strip.background = element_blank(), 
          strip.text = element_text(hjust=0, face="bold"),
          aspect.ratio = .7,
          panel.grid = element_blank(), 
          legend.position=c(.955,.37), 
          legend.background = element_blank(), 
          legend.key.height=unit(.9,"line")
    ) +
    scale_y_continuous(breaks = seq(-4, 4, 2), expand = expansion(mult=.03)) +
    scale_x_continuous(expand = expansion(mult=.03)) +
    guides(fill="none", colour="none", lty="none") +
    labs(x = "Scaled elevation", 
         y = "Log-odds ratio", colour = "", fill= "") 

## save ----
p_occ_all <- egg::ggarrange(occ_p1, occ_p2, occ_p3)
ggsave("figures/hyperparameter_probability_and_logOR.png", plot=p_occ_all, 
       width=249*.8, height=121*1.6, 
       units="mm", dpi=400)


# species level plots ----
newdat_sp <- expand.grid(range_pos = seq(-1, 1, len=30), 
                         species = unique(fit$data$species), 
                         habitat = c(1, .6), 
                         # medium_dependency = c(0, 1), 
                         time = 1, 
                         obsvr_point = "JS") %>%
    left_join(., unique(dplyr::select(analysis_df, species, 
                                      ele_midpoint=midpoint_both, 
                                      scale_both))) %>%
    left_join(., unique(dplyr::select(fit$data, species, species_eltontraits, 
                                      medium_dependency))) %>%
    mutate(range_pos2 = range_pos^2, 
           range_pos_upr = ifelse(range_pos < 0, 0, range_pos),
           range_size_sc = scale(scale_both),
           ele_midpoint_std = (ele_midpoint - mean(analysis_df$midpoint_both))/sd(analysis_df$midpoint_both), 
           habitat_std = (habitat - mean(analysis_df$habitat))/sd(analysis_df$habitat), 
           ele_pt = scale_both * range_pos + ele_midpoint, 
           ele_pt_std = (ele_pt - mean(analysis_df$ele_pt))/sd(analysis_df$ele_pt))

preds_sp <- fitted_flocker(fit, type = "occupancy", 
                           new_data = newdat_sp, 
                           re_formula = ~ (1 + habitat_std + range_pos + range_pos2|species) + 
                               (1 + habitat_std|gr(species_eltontraits, cov=A)), 
                           summarise = F, ndraws = 500, response = F)

## calculate log odds ----
logOdds_sp <- cbind(newdat_sp, preds_sp) %>%
    group_by(ele_midpoint, range_pos, ele_pt, medium_dependency, species) %>%
    summarise_at(vars(paste0("iter_", 1:500)), .funs = function(x) x[2]-x[1])

logOdds2_sp <- logOdds_sp %>%
    ungroup %>%
    dplyr::select(iter_1:iter_500) %>%
    as.matrix(.)

logOdds3_sp <- logOdds_sp %>%
    ungroup %>%
    mutate(mean = rowMeans(logOdds2_sp), 
           lwr = matrixStats::rowQuantiles(logOdds2_sp, probs = .05),
           upr = matrixStats::rowQuantiles(logOdds2_sp, probs = .95)) %>%
    dplyr::select(-(iter_1:iter_500)) %>%
    mutate(label = case_when(medium_dependency == dep_lvls["high"] ~ "(c) Odds ratio (high dependency)", 
                             medium_dependency == dep_lvls["medium"] ~ "(d) Odds ratio (medium dependency)"))

## summarise occupancy ----
median_occ_sp <- cbind(newdat_sp, preds_sp) %>%
    group_by(ele_midpoint, range_pos, ele_pt, medium_dependency, species) %>%
    summarise_at(vars(paste0("iter_", 1:500)), .funs = mean)

median_occ_mat <- median_occ_sp %>%
    ungroup %>%
    dplyr::select(iter_1:iter_500) %>%
    as.matrix(.) %>%
    boot::inv.logit(.)

median_occ_sp2 <- median_occ_sp %>%
    ungroup %>%
    mutate(mean = rowMeans(median_occ_mat), 
           lwr = matrixStats::rowQuantiles(median_occ_mat, probs = .05),
           upr = matrixStats::rowQuantiles(median_occ_mat, probs = .95)) %>%
    dplyr::select(-(iter_1:iter_500)) %>%
    mutate(label = case_when(medium_dependency == dep_lvls["high"] ~ "(a) Occupancy (high dependency)", 
                             medium_dependency == dep_lvls["medium"] ~ "(b) Occupancy (medium dependency)"))

## plot ----
p1 <- ggplot(median_occ_sp2, aes(ele_pt, mean, group=species, colour = ele_midpoint, 
                           lty = factor(medium_dependency))) + 
    facet_wrap(~label) +
    geom_line(alpha=.8) +  
    scale_colour_viridis_c() +
    scale_linetype_manual(values = c(2,1)) +
    coord_cartesian(xlim=c(880, 3800)) +
    labs(x = "Elevation (m)", y = "P(Occupancy)") +
    theme_bw() +
    theme(axis.text = element_text(colour="black"), 
          strip.text = element_text(face="bold", hjust=0), 
          strip.background = element_blank(),
          panel.grid = element_blank()#,
          #aspect.ratio = .3
    ) + 
    guides(lty="none", col="none")

p2 <- ggplot(logOdds3_sp, aes(ele_pt, mean, group=species, col = ele_midpoint, 
                              lty = factor(medium_dependency))) + 
    facet_wrap(~label) +
    geom_line(alpha=.8) +  
    geom_hline(yintercept = 0, lty="longdash") +
    scale_colour_viridis_c() +
    scale_linetype_manual(values = c(2,1)) +
    coord_cartesian(xlim=c(880, 3800)) +
    labs(x = "Elevation (m)", y = "Log-odds ratio") +
    theme_bw() +
    theme(axis.text = element_text(colour="black"), 
          strip.text = element_text(face="bold", hjust=0), 
          strip.background = element_blank(),
          panel.grid = element_blank()#,
          #aspect.ratio = .3
    ) + 
    guides(lty="none", col="none")

## save ----
plot_species_both <-egg::ggarrange(p1  + 
                                       theme(axis.text.x = element_blank(), 
                                             axis.title.x=element_blank()), p2, 
                                   ncol=1)

ggsave("figures/species_probablity_logOR.png", plot = plot_species_both,
       units="mm", width=190, height=125, dpi=400)


# species occupancy at midpoint ----
## high dependency ----
p1 <- logOdds3_sp %>%
    filter(medium_dependency == 1) %>%
    filter(range_pos == sort(unique(range_pos))[15]) %>%
    arrange(mean) %>%
    mutate(species = gsub("_", " ", species), 
           species = factor(species, levels=species)) %>%
    filter(species %in% unique(species)[1:50]) %>%
    ggplot(aes(species, mean)) + 
    geom_linerange(aes(ymin=lwr, ymax=upr), 
                   position = position_dodge(width=.8)) +
    geom_point(, position = position_dodge(width=.8)) +
    geom_hline(yintercept=0, linetype="longdash") +
    scale_shape_manual(values=c(16, 21)) +
    scale_fill_manual(values=c("black", "white")) +
    labs(x="", y="Log-odds ratio", shape = "Model", fill="Model") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position=c(.95,.2), 
          legend.background = element_rect(colour="black"), 
          plot.margin = unit(c(1, 1, 1, 10), "mm")) 

p2 <- logOdds3_sp %>%
    filter(medium_dependency == 1) %>%
    filter(range_pos == sort(unique(range_pos))[15]) %>%
    arrange(mean) %>%
    mutate(species = gsub("_", " ", species), 
           species = factor(species, levels=species)) %>%
    filter(species %in% unique(species)[51:100]) %>%
    ggplot(aes(species, mean)) + 
    geom_linerange(aes(ymin=lwr, ymax=upr), 
                   position = position_dodge(width=.8)) +
    geom_point(position = position_dodge(width=.8)) +
    geom_hline(yintercept=0, linetype="longdash") +
    scale_shape_manual(values=c(16, 21)) +
    scale_fill_manual(values=c("black", "white")) +
    labs(x="", y="Log-odds ratio", shape = "Model", fill="Model") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position=c(.95,.2), 
          legend.background = element_rect(colour="black"), 
          plot.margin = unit(c(1, 1, 1, 10), "mm")) 

## high dep save ----
plot_OR_range_centre <- egg::ggarrange(p2, p1, ncol=1)
ggsave("figures/OR_range_centre_high_dep.png", plot_OR_range_centre, 
       height=210, width=293, units="mm")

## medium dependency ----
p1 <- logOdds3_sp %>%
    filter(medium_dependency == -1) %>%
    filter(range_pos == sort(unique(range_pos))[15]) %>%
    arrange(mean) %>%
    mutate(species = gsub("_", " ", species), 
           species = factor(species, levels=species)) %>%
    filter(species %in% unique(species)[1:71]) %>%
    ggplot(aes(species, mean)) + 
    geom_linerange(aes(ymin=lwr, ymax=upr), 
                   position = position_dodge(width=.8)) +
    geom_point(, position = position_dodge(width=.8)) +
    geom_hline(yintercept=0, linetype="longdash") +
    scale_shape_manual(values=c(16, 21)) +
    scale_fill_manual(values=c("black", "white")) +
    labs(x="", y="Log-odds ratio", shape = "Model", fill="Model") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position=c(.95,.2), 
          legend.background = element_rect(colour="black"), 
          plot.margin = unit(c(1, 1, 1, 10), "mm")) 

p2 <- logOdds3_sp %>%
    filter(medium_dependency == -1) %>%
    filter(range_pos == sort(unique(range_pos))[15]) %>%
    arrange(mean) %>%
    mutate(species = gsub("_", " ", species), 
           species = factor(species, levels=species)) %>%
    filter(species %in% unique(species)[72:143]) %>%
    ggplot(aes(species, mean)) + 
    geom_linerange(aes(ymin=lwr, ymax=upr), 
                   position = position_dodge(width=.8)) +
    geom_point(, position = position_dodge(width=.8)) +
    geom_hline(yintercept=0, linetype="longdash") +
    scale_shape_manual(values=c(16, 21)) +
    scale_fill_manual(values=c("black", "white")) +
    labs(x="", y="Log-odds ratio", shape = "Model", fill="Model") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          panel.grid=element_blank(), 
          axis.text = element_text(colour="black"), 
          legend.position=c(.95,.2), 
          legend.background = element_rect(colour="black"), 
          plot.margin = unit(c(1, 1, 1, 10), "mm")) 

## medium dep save ----
plot_OR_range_centre <- egg::ggarrange(p2, p1, ncol=1)
ggsave("figures/OR_range_centre_med_dep.png", plot_OR_range_centre, 
       height=210, width=293, units="mm")


# fixed effects ----
feffs <- fixef(fit, summary=F) %>%
    as_tibble %>%
    slice(1:1000)

high_dep_effs <- feffs %>%
    as_tibble() %>%
    mutate(inter_eff = `occ_ele_midpoint_std:habitat_std`, 
           intra_eff = `occ_range_pos:habitat_std`, 
           intra_eff_upper = `occ_range_pos:habitat_std` + 
               `occ_habitat_std:range_pos_upr`) %>%
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
    mutate(inter_eff = `occ_ele_midpoint_std:habitat_std` + 
               `occ_ele_midpoint_std:habitat_std:medium_dependency`, 
           intra_eff = `occ_range_pos:habitat_std` + 
               `occ_range_pos:habitat_std:medium_dependency`, 
           intra_eff_upper = intra_eff + `occ_habitat_std:range_pos_upr` + 
               `occ_habitat_std:medium_dependency:range_pos_upr`) %>%
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
          aspect.ratio = .7,
          panel.grid = element_blank(), 
          legend.position=c(.86,.85), 
          legend.background = element_rect(colour="black"), 
          legend.key.height=unit(.9,"line"), 
          plot.margin = unit(c(0.2, 2, 0.2, 0.2), units = "cm")
    ) +
    #guides(fill="none") +
    geom_hline(yintercept = 0, lty="longdash") +
    labs(y = "", x = "", fill = "Dependency") +
    geom_text(aes(label=pd, y = upr), position = position_dodge(width=.5), 
              hjust=-0.2, vjust=.3, size=3)


## calculate lambda ----
sigmas <- VarCorr(fit, summary=F)
lambdas <- (sigmas$species_eltontraits$sd/(sigmas$species$sd[,3:4] + 
                                               sigmas$species_eltontraits$sd))[1:1000,] %>%
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
    filter(par == "occ_habitat_std") %>%
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

## save ----
effs_both <- egg::ggarrange(effs_p1, effs_p2, ncol=1, heights=c(1, .1))
ggsave("figures/parameter_plot.png", plot = effs_both)
