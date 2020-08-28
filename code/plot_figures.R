## Plot figures (that don't depend on model output)

# depends:
#       df_analysis
#       df_pointInfo
# note: need to also save zf_subset.. 

library(sf); library(dplyr); library(ggplot2); library(egg)


## number parameters 
# *not* hyperpars
sum(stan_data$n_species*6, stan_data$n_cluster_species, stan_data$n_site_species)
sum(stan_data$n_species*7, stan_data$n_cluster_species, stan_data$n_site_species)


## Point elevations & configuration ----
points_sp <- st_as_sf(df_ptInfo, coords=c("lon", "lat"), crs="epsg:4326") %>%
    st_transform(., crs="epsg:3117") %>%
    bind_cols(., as_tibble(st_coordinates(.))) %>%
    rename(lon_trans = `X`, lat_trans=`Y`) %>%
    group_by(rotation) %>%
    mutate(lon_trans = lon_trans - min(lon_trans), 
           lat_trans = lat_trans - min(lat_trans)) %>%
    group_by(rotation) %>%
    mutate(midpoint = median(ele_jaxa)) %>%
    ungroup %>%
    arrange(desc(midpoint)) %>%
    mutate(rot = gsub("rot_", "", rotation), 
           rot=factor(rot, levels = unique(rot)))

st_coordinates(points_sp) %>%
    apply(., 2, function(x) diff(range(x))/1000)

max_dist <- st_coordinates(points_sp) %>%
    dist() %>% 
    max
max_dist/1000

#* plots ----

breaks_fun <- function(x) {y <- round(max(x), 0); c(0, y - y %% 10^(nchar(y)-1))}

p1 <- ggplot(points_sp, aes(lon_trans, lat_trans, fill=ele_jaxa)) + 
    geom_point(pch=21) + 
    scale_fill_gradientn(colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", 
                                     "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", 
                                     "#D73027")) +
        scale_x_continuous(breaks=breaks_fun, expand=expansion(add=300)) +
        scale_y_continuous(breaks=breaks_fun, expand=expansion(add=300)) +
    facet_wrap(~rot, scales="free", ncol=5) +
    theme(strip.background = element_blank(), 
          strip.text = element_text(face="bold", hjust=0), 
          axis.text.y = element_text(colour="black", angle=45), 
          panel.grid.minor = element_blank(), 
          axis.title = element_blank(),
          axis.text.x = element_text(colour="black", angle=45, hjust=1, vjust=1.1)) +
    guides(fill=F) +
    labs(x="", y="")

p2 <- ggplot(points_sp, aes(rot, ele_jaxa, fill=ele_jaxa)) + 
    geom_jitter(width=.05, height=0, pch=21) + 
    scale_fill_gradientn(colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                                     "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43",
                                     "#D73027")) +
    scale_y_continuous(breaks=seq(0, 4000, 1000), limits=c(0, 4000), expand=c(0,0)) +
    # scale_x_discrete(expand=c(0,0.3)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5), 
          aspect.ratio=.8, 
          axis.text=element_text(colour="black"), 
          panel.grid.major = element_line(linetype = 1, colour="grey95")) +
    labs(y = "Elevation (m)", x="Site") +
    guides(fill=F)

ggsave("figures/points_by_rotation.png", p1, width=150, height=200, units="mm")
ggsave("figures/elevation_by_site.png", p2)

## Pairwise distances ----
dist_as_df <- function(df_i) {
    id_pt <- df_i$point
    id_cl <- df_i$cluster
    names(id_cl) <- id_pt
    d <- dist(st_coordinates(df_i)) %>%
        as.numeric
    pt_combn <- combn(id_pt, m=2) 
    tibble(x = pt_combn[1,], y = pt_combn[2,]) %>%
        mutate_all(as.character) %>%
        mutate(x_cl = recode(x, !!!id_cl), 
               y_cl = recode(y, !!!id_cl), 
               within_cl = x_cl == y_cl,
               d = d)    
}

# get distances
dists <- points_sp %>%
    split(., .$rotation) %>%
    lapply(., dist_as_df) %>%
    bind_rows(., .id="rot")

plot_cl_dist_1 <- ggplot(dists, aes(d, fill=within_cl)) + 
    geom_histogram(position="identity", alpha=.5, binwidth=100, boundary=0) +
    labs(x="Distance (m)", y="Frequency", fill="Within cluster") +
    scale_x_continuous(expand=c(0.01, 0)) +
    scale_y_continuous(expand=c(0.01, 0)) +
    guides(fill=F) +
    theme_classic() +
    theme(aspect.ratio = .7, 
          axis.title.x = element_blank(),
          axis.line.y = element_line(size=1))

dist_lims <- c(0, max(dists$d))
plot_cl_dist_2 <- ggplot(dists, aes(d, col=within_cl)) + 
    geom_vline(xintercept=c(200, 500), lty=2, col="grey40") +
    stat_ecdf() +
    # geom_histogram(position="identity", alpha=.5, binwidth=100, boundary=0) +
    labs(x="Distance (m)", y="Cumulative proportion", fill="Within cluster") +
    scale_x_continuous(expand=c(0.01, 0), limits=dist_lims) +
    scale_y_continuous(expand=c(0.01, 0)) +
    guides(col=F) +
    theme_classic() +
    theme(aspect.ratio = .7, 
          axis.line.y = element_line(size=1),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey95"))

mean_dists <- dists %>% group_by(within_cl) %>%
    summarise(m = median(d), 
              lwr = quantile(d, .95), 
              upr = quantile(d, .05), 
              min = min(d), 
              max = max(d))

plot_cl_dist_3 <- ggplot(mean_dists, aes(m, y=within_cl, xmin=lwr, xmax=upr, 
                                         col=within_cl)) + 
    geom_point() +
    geom_errorbarh(height=0, size=1) +
    geom_errorbarh(aes(xmin=min, xmax=max), height=0, alpha=.5, size=1) +
    scale_x_continuous(limits=dist_lims) +
    labs(x="Distance (m)") +
    guides(col=F) +
    theme_classic() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank())

plot_cl_dist_all <- ggarrange(plot_cl_dist_1, plot_cl_dist_2, plot_cl_dist_3, 
                              ncol=1, heights=c(1, 1, .1))

ggsave("figures/plot_clusterdists.png", plot_cl_dist_all)

## Species:point detections ----
points_by_species <- analysis_df %>% 
    filter(Q == 1) %>%
    group_by(species) %>%
    summarise(n_Q = n()) 

# 52% of species detected on fewer than 4 points
points_by_species %>% summarise(sum(n_Q < 4)/n())
points_by_species %>% summarise(sum(n_Q < 3)/n())
points_by_species %>% summarise(sum(n_Q < 2)/n())
points_by_species %>% summarise(sum(n_Q >= 10)/n())

plot_points_by_species1 <- ggplot(points_by_species, aes(n_Q)) + stat_ecdf() +
    scale_x_continuous(breaks=c(2, 4, 8, 16, 32, 64)) +
    theme_classic() +
    labs(x = "Number of point-detections per species", 
         y = "Cumulative proportion")

plot_points_by_species1 <-ggplot(points_by_species, aes(n_Q, y=cumsum(..count..)/sum(..count..))) + 
    geom_histogram(boundary=0.5, binwidth=1, fill="grey90", colour="black") +
    scale_x_continuous(breaks=c(2, 4, 8, 16, 32, 64, 128)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(colour="black"),
          axis.line.y = element_line(colour="black", size=1),
          panel.grid.major.y = element_line(colour="grey95"), 
          aspect.ratio = .5) +
    labs(x = "Number of point-detections per species", 
         y = "Cumulative proportion")

plot_points_by_species2 <-ggplot(points_by_species, aes(n_Q, y=..count../sum(..count..))) + 
    geom_histogram(boundary=0.5, binwidth=1, fill="grey90", colour="black") +
    scale_x_continuous(breaks=c(2, 4, 8, 16, 32, 64, 128)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0, .3, .1)) +
    theme_classic() +
    theme(axis.text = element_text(colour="black"), 
          axis.line.y = element_line(colour="black", size=1),
          panel.grid.major.y = element_line(colour="grey95"), 
          aspect.ratio = .5) +
    labs(x = "Number of point-detections per species", 
         y = "Proportion") +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank())

plot_points_by_species_both <- egg::ggarrange(plot_points_by_species2, 
                                              plot_points_by_species1, 
                                              ncol=1)
ggsave("figures/density_points_by_species.png", 
       plot = plot_points_by_species_both)


## Elevation histogram ----
eleRange_lims <- analysis_df$ele_sc %>% range
eleHist_1 <- ggplot(analysis_df, aes(ele_sc, fill=factor(Q), group=factor(Q))) + 
    geom_histogram(binwidth=.2, boundary=0, col="black") +
    labs(x="Scaled elevation", y="Frequency") +
    scale_y_continuous(expand = expansion(0,0), breaks=seq(0,2000, 500)) +
    scale_x_continuous(expand = expansion(0), breaks=c(-5, -2, -1, 0, 1, 2, 5), 
                       limits = eleRange_lims) +
    scale_fill_manual(values=c("grey90", "grey50")) +
    theme(panel.background = element_blank(), 
          #aspect.ratio=.7,
          axis.line.y = element_line(),
          panel.grid = element_blank(),
          axis.text=element_text(colour="black"), 
          axis.title.x = element_blank(),
          axis.text.x=element_blank()) + 
    guides(fill=FALSE) #+
    # geom_rangeframe(data=data.frame(x=c(-5,5), y=c(0,2000)), aes(x, y), inherit.aes = F)

props <- analysis_df %>% mutate(ele_f = cut(ele_sc, seq(-50, 50, .2), 
                                            labels = seq(-49.9, 49.9, .2))) %>%
    group_by(Q, ele_f) %>%
    summarise(N = n()) %>%
    group_by(ele_f) %>% 
    mutate(N_tot = sum(N)) %>%
    summarise(N_rel = N[2]/N_tot) %>%
    mutate(ele_f = as.numeric(as.character(ele_f)))


eleHist_2 <- ggplot(props, aes(ele_f, N_rel)) + 
    # geom_col(col="black", width = .1,) +
    stat_summary(fun=sum,geom="bar", colour="black", fill="grey50", width=.2) +
    scale_y_continuous(expand = expansion(0)) +
    scale_x_continuous(expand = expansion(0), breaks=c(-5, -2, -1, 0, 1, 2, 5), 
                       limits = eleRange_lims) +
    theme(panel.background = element_blank(), 
          axis.line.x = element_line(), 
          axis.line.y = element_line(),
          panel.grid = element_blank(),
          axis.text=element_text(colour="black")) +
    labs(x="Scaled elevation", y="Proportion") #+
    # geom_rangeframe(data=data.frame(x=c(-5,5), y=c(0,.2)), aes(x, y))

p_ele_both <- egg::ggarrange(eleHist_1, eleHist_2, heights=c(1,.3))
ggsave("figures/elevational_histograms.png", p_ele_both, width=200, height=150, units="mm")


## Ele limit correspondence----
rangeMar_p1 <- ggplot(ele_both, aes(lwr_Quinones, lwr_McMullan)) + 
    geom_jitter(height=50, width=50) + 
    coord_fixed() +
    geom_abline() +
    theme_classic() +
    theme(axis.text = element_text(colour="black"), 
          axis.line.y = element_line(colour="black", size = 1)) +
    labs(x="Lower margin (Quiñones)", y="Lower margin (McMullan)")

rangeMar_p2 <- ggplot(ele_both, aes(upr_Quinones, upr_McMullan)) + 
    geom_jitter(height=50, width=50) + 
    coord_equal() +
    geom_abline() +
    theme_classic() +
    scale_x_continuous(expand=c(0,50)) +
    scale_y_continuous(expand=c(0,50)) +
    theme(axis.text = element_text(colour="black"), 
          axis.line.y = element_line(colour="black", size = 1)) +
    labs(x="Upper margin (Quiñones)", y="Upper margin (McMullan)")

rangeMar_both <- ggarrange(rangeMar_p1, rangeMar_p2, ncol=2)
ggsave("figures/rangeMargin_comparison.png", plot=rangeMar_both)

cor(ele_both$lwr_McMullan, ele_both$lwr_Quinones, use = "pairwise.complete.obs")
cor(ele_both$upr_McMullan, ele_both$upr_Quinones, use = "pairwise.complete.obs")

## time of day ----
c(round(min(zf_subset$time)/60, 0), round(((min(zf_subset$time)/60) %% 1)*60, 0))

c(round(max(zf_subset$time)/60, 0), round(((max(zf_subset$time)/60) %% 1)*60, 0))

visit_times <- zf_subset %>% select(point, visit, time) %>% unique 

ggplot(visit_times, aes(time/60)) + 
    geom_histogram(boundary=0, binwidth = .25, fill="grey90", col="black") +
    geom_vline(xintercept=c(5.5, 12), lty=2) +
    scale_x_continuous(breaks=seq(0, 14, 1)) +
    theme_classic() +
    theme(aspect.ratio = .6) +
    labs(x="Time of day (hours)", 
         y = "Frequency")
ggsave("figures/visit_times.png")
# 15 counts, 2% run later in day than 12
visit_times %>% summarise(sum(time/60 > 12))
visit_times %>% summarise(sum(time/60 > 12)/n())

5.45 %% 1
min(zf_subset$time)/60
5.45 %% 1
