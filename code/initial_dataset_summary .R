## Create initial summary of EC dataset
## Need to understand: 
#   - Number of observations per species
#   - Number obs per species and habitat type
#   - Broad elevational associations (i.e. when scaled how many of them neatly 
#   follow unimodal distribution etc.)
#   - Variation in habitat metrics by elevation (not done yet)

# housekeeping ----
library(dplyr); library(boot)

df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds")
unique_pt <- df_bird$Point %>% unique
unique_ptVis <- df_bird %>%
    group_by(Point) %>%
    summarise(nVis = length(unique(Visit)))

## forest cover and landscape metrics ----
forest_cover <- readRDS("data/hansen_500m_buffer.rds") %>%
    mutate(tc_2000 = treecover2000 > 50, 
           tc_2018 = as.numeric(tc_2000 - (lossyear != 0) == 1)) %>%
    group_by(point_id) %>%
    summarise(amt_tc = sum(tc_2018), pct_tc = amt_tc/n())

lscape_metrics <- readRDS("data_landscapeMetrics_draft_v1.rds") %>%
    group_by(Point) %>% slice(1)

## elevational ranges (Quinones)
ele_range <- read.csv("../Colombia/data/bird data_Jan&Jun2019/elevational_ranges_Quinones.csv") %>%
    mutate(Species = paste0(Genus, "_", Species)) %>%
    mutate(Lower = ifelse(is.na(Lower), 0, Lower))

## point info 
df_ptInfo <- read.csv("data/CO_sampling_points_metafile_ver2.csv") %>%
    as_tibble() %>%
    left_join(., forest_cover) %>%
    mutate(hab = substr(point_id, 3, 3)) %>%
    rename(Point = point_id) %>%
    left_join(., lscape_metrics) %>%
    filter(Point %in% unique_pt) %>%
    select(Point, lat, long, ele=ALOSelev, amt_tc, pct_tc, Habitat, 
           pct_nat_hab = Prop_natural_habitat)


## get full dataframe and zero-filled version
df_full <- left_join(df_bird, df_ptInfo) %>%
    left_join(., ele_range) %>%
    #filter(hab == "P") %>%
    group_by(Point, Species) %>%
    mutate(Z = 1) %>%
    slice(1) %>%
    select(Site, Point, Species, Z, lat:Z)

## Number of points by species ----
nObs_bySpecies <- df_full %>% 
    group_by(Species) %>% 
    summarise(N=n()) %>%
    arrange(N)

nObs_bySp_forest <- df_full %>% 
    filter(Habitat == "Forest") %>%
    group_by(Species) %>% 
    summarise(N=n()) %>%
    arrange(N)

nObs_bySp_pasture <- df_full %>% 
    filter(Habitat == "Pasture") %>%
    group_by(Species) %>% 
    summarise(N=n()) %>%
    arrange(N)

lab_all <- paste0("#Sp(N>10) = ", sum(nObs_bySpecies$N >=10), "\n", 
                  "#Sp(N>20) = ", sum(nObs_bySpecies$N >=20), "\n",  
                  "#Sp(N>40) = ", sum(nObs_bySpecies$N >=40), "\n", 
                  "#Sp(N>80) = ", sum(nObs_bySpecies$N >=80))

lab_F <- paste0("#Sp(N>10) = ", sum(nObs_bySp_forest$N >=10), "\n", 
                "#Sp(N>20) = ", sum(nObs_bySp_forest$N >=20), "\n",  
                "#Sp(N>40) = ", sum(nObs_bySp_forest$N >=40), "\n", 
                "#Sp(N>80) = ", sum(nObs_bySp_forest$N >=80))

lab_P <- paste0("#Sp(N>10) = ", sum(nObs_bySp_pasture$N >=10), "\n", 
                "#Sp(N>20) = ", sum(nObs_bySp_pasture$N >=20), "\n",  
                "#Sp(N>40) = ", sum(nObs_bySp_pasture$N >=40), "\n", 
                "#Sp(N>80) = ", sum(nObs_bySp_pasture$N >=80))


p1 <- ggplot(nObs_bySpecies, aes(N)) + stat_bin(binwidth = 1, aes(y=cumsum(..count..)), geom="step")+ 
    scale_x_continuous(breaks=c(0, 5, 10, 20, 40, 80)) +
    scale_y_continuous(breaks=seq(0, 550, 50), 
                       labels = paste0(seq(0, 550, 50), " (", round(seq(0, 550, 50)/nrow(nObs_bySpecies), 2)*100, "%)")) +
    theme(panel.grid.minor = element_blank()) +
    geom_label(x=Inf, y=-Inf, label=lab_all, hjust=1.1, vjust=-.1)


p2 <- ggplot(nObs_bySp_forest, aes(N)) + 
    stat_bin(binwidth = 1, aes(y=cumsum(..count..)), geom="step")+ 
    scale_x_continuous(breaks=c(0, 5, 10, 20, 40, 80)) +
    scale_y_continuous(breaks=seq(0, 550, 50), 
                       labels = paste0(seq(0, 550, 50), " (", round(seq(0, 550, 50)/nrow(nObs_bySp_forest), 2)*100, "%)")) +
    theme(panel.grid.minor = element_blank()) +
    geom_label(x=Inf, y=-Inf, label=lab_F, hjust=1.1, vjust=-.1) +
    labs(title="Forest-points only")

p3 <- ggplot(nObs_bySp_pasture, aes(N)) + 
    stat_bin(binwidth = 1, aes(y=cumsum(..count..)), geom="step")+ 
    scale_x_continuous(breaks=c(0, 5, 10, 20, 40, 80)) +
    scale_y_continuous(breaks=seq(0, 550, 50), 
                       labels = paste0(seq(0, 550, 50), " (", round(seq(0, 550, 50)/nrow(nObs_bySp_forest), 2)*100, "%)")) +
    theme(panel.grid.minor = element_blank()) +
    geom_label(x=Inf, y=-Inf, label=lab_P, hjust=1.1, vjust=-.1) +
    labs(title="Pasture-points only")

p_all <- egg::ggarrange(p1, p2, p3, ncol=1)
ggsave("figures/nPoints_bySpecies.png", plot=p_all, width=6, height=12)


## obs by elevation
species_subset <- df_full %>%
    group_by(Species) %>%
    summarise(N=n()) %>%
    filter(N >= 10)

df_subset <- df_full %>% 
    filter(Species %in% species_subset$Species) %>% 
    dplyr::select(Point, Species, Z)

zf_subset <- expand.grid(Species = species_subset$Species, Point = unique(df_full$Point)) %>%
    mutate(Z = 0) %>%
    bind_rows(., df_subset) %>%
    group_by(Point, Species) %>%
    summarise(Z = max(Z)) %>%
    left_join(., df_ptInfo) %>%
    left_join(., ele_range) %>%
    mutate(ele_sc = (ele- (Upper+Lower)/2)/(Upper-Lower)) %>%
    filter(!grepl("D$", Point)) %>%
    filter(!is.na(ele_sc)) %>%
    mutate(Site_sp= paste0(substr(Point, 1,2), "_", Species))

obs_byEle_full <- zf_subset %>% 
    mutate(ele_f = cut(ele, seq(0, 5000, 100), labels = seq(50, 4950, 100)), 
           ele_f = as.numeric(as.character(ele_f))) 

obs_byEle <- obs_byEle_full %>%
    group_by(Species, ele_f) %>%
    summarise(sumZ = sum(Z), N = n(), ele_sc = unique((ele_f- (Upper+Lower)/2)/(Upper-Lower)))


sp_list <- obs_byEle$Species %>% unique
pred_list <- list()
library(MASS)
for(i in 1:length(sp_list)) {
    sp_i <- zf_subset %>% filter(Species == sp_list[i])
    fit_i <- glm(Z ~ poly(ele_sc, 2), sp_i, family = "binomial")
    # fit_i <- glm(cbind(sumZ, N) ~ poly(ele_sc, 2), sp_i, family = "binomial")
    pred_list[[i]] <- with(sp_i, tibble(Species = unique(Species), 
                                ele_sc = seq(min(ele_sc), max(ele_sc), len=50))) %>%
        mutate(p = predict(fit_i, newdata=.))
}
pred_all <- bind_rows(pred_list)

sp_samp <- sample(species_subset$Species, 22)
shaded_area <- obs_byEle %>% group_by(Species) %>% summarise(ymax = max(sumZ/N))
ggplot(obs_byEle %>% filter(Species %in% sp_samp), aes(ele_sc, sumZ/N)) + 
    annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, alpha = .3) +
    geom_point() + 
    facet_wrap(~Species, scales="free_y", ncol=3) +
    geom_vline(xintercept=0, lty=3) +
    geom_line(data=pred_all %>% filter(Species %in% sp_samp), aes(y=boot::inv.logit(p)), col="red") +
    theme_bw() +
    theme(strip.text = element_text(hjust=0, face="bold"), 
          strip.background = element_blank(), 
          axis.text = element_text(colour="black"), 
          panel.grid= element_blank()) +
    labs(x = "Standardised elevation", y="Occupancy")
ggsave("figures/elevational_distributions.png", height=200, width=150, units="mm")


ele_summ <- pred_all %>% group_by(Species) %>%
    summarise(ele_maxP = ele_sc[which.max(p)]) %>%
    left_join(., ele_range) %>%
    mutate(mid = (Lower+Upper)/2, ele_unsc = ele_maxP * (Upper-Lower) + (Lower+Upper)/2)

ggplot(ele_summ, aes(mid, ele_unsc)) + geom_point() + geom_abline() + coord_equal() +
    labs(x="Elevational range mid-point (Quiñones)", y="Elevation with max(abundance)", 
         caption = "Note: Estimated elevational midpoints for species with N\u226510")
ggsave("figures/elevational_midpoints.png")


# get example plots
# forest_cover %>%
#     # as.data.frame %>% 
#     filter(id %in% sample(unique(id), 20)) %>% 
#     mutate(tc = cut(treecover2000, c(-.1, 25, 50, 75, 100.1), ordered=T)) %>%
#     ggplot(aes(x, y, fill=tc_2018)) + geom_tile() + 
#     facet_wrap(~point_id, scales="free") +
#     theme_bw() +
#     theme(aspect.ratio=1, 
#           strip.text=element_text(hjust=0, face="bold"), 
#           strip.background = element_blank(), 
#           panel.border = element_blank(), 
#           panel.grid=element_blank(), 
#           axis.text = element_blank(), 
#           axis.ticks=element_blank()) +
#     labs(x="", y="")
# ggsave("figures/example_tc2000.png")