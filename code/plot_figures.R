# ## Plot figures that don't depend on model output
# 
# # depends:
# #       df_analysis
# #       df_pointInfo
# # note: need to also save zf_subset.. 
# #
# # 
# analysis_df_full <- left_join(det_vis, df_ptInfo) %>%
#   left_join(., ele_both) %>%
#   left_join(., obsvr) %>%
#   ungroup %>%
#   #filter(species %in% sample(unique(species), 100)) %>%
#   mutate(#ele_sc = round(ele-midpoint_both, 0)/2000, 
#     ele_sc = range_scaling(ele, midpoint_both, scale_both), 
#     habitat_sc = (habitat - mean(habitat))/sd(habitat),
#     # half = ifelse(ele_sc < 0, 0, 1), 
#     # dist = ele_sc^2, 
#     # dist_upr = ifelse(ele_sc < 0, 0, dist), 
#     # hab_dist = habitat_sc * dist, 
#     # hab_dist_upr = ifelse(ele_sc < 0, 0, hab_dist),
#     id_sp = as.numeric(as.factor(species)),
#     id_pt = as.numeric(as.factor(point)),
#     id_site_sp = paste0(rotation, species),
#     id_site_sp_num = as.numeric(as.factor(id_site_sp)), 
#     id_cl_sp = paste0(cluster, species),
#     id_cl_sp_num = as.numeric(as.factor(id_cl_sp)), 
#     Q = as.numeric(rowSums(select(., d1:d4), na.rm = T)>0), 
#     n_visit = rowSums((select(., d1:d4) != -99))) %>%
#   # tidy up a bit (drop redundant cols) 
#   select(-donegan, -donegan2, -ebird2, -hbw, -ebird, -pulido, -eltontraits, 
#          -english, -genus, -is_na) 
# ### 
library(sf); library(dplyr); library(ggplot2); library(egg)

analysis_df <- readRDS("output/analysis_df_birdlife_fDep_High.rds")

## number parameters 
# *not* hyperpars
sum(stan_data$n_species*6, stan_data$n_cluster_species, stan_data$n_site_species)
sum(stan_data$n_species*7, stan_data$n_cluster_species, stan_data$n_site_species)

# Supp table 1 ----
for_write <- analysis_df %>%
  filter(Q==1) %>%
  mutate(lwr_both = midpoint_both - scale_both, upr_both = midpoint_both + scale_both) %>%
  group_by(species_hbw) %>%
  #filter(ele %in% c(min(ele), max(ele))) %>%
  group_by(species_hbw, lwr_both, upr_both) %>% 
  summarise(species = unique(gsub("_", " ", species_hbw)), 
            max_ele = max(ele), min_ele = min(ele), 
            max_range_pos = round(max(ele_sc), 2), 
            min_range_pos = round(min(ele_sc),2),
            max_m_from_upper = unique(max_ele - upr_both), 
            min_m_from_upper = unique(lwr_both - min_ele),
            beyond_upper = unique(max_ele > upr_both), 
            beyond_lower = unique(min_ele < lwr_both), 
            either = beyond_upper|beyond_lower, 
            maxDiff = max(c(abs(max_range_pos), abs(min_range_pos)))) %>%
  arrange(desc(maxDiff))

write.csv(for_write, "output/species_summary_table.csv")

# number of points
length(unique(analysis_df$point))

# range of elevations 
range(analysis_df$ele)

## Landscape correlations ----
tc_metrics <- readRDS("data/forest_cover_metrics_500m.rds") %>%
  left_join(., readRDS("data/point_ele_EasternCordillera.rds")) %>%
  filter(!grepl("^CC", point) & ele_jaxa >= 1000) %>%
  select(point, amt_tc, pct_tc, edge, contig, pct_core_area, dist_from_edge)

# correlations
cor(tc_metrics[,-1], method="spearman")

tc_metrics %>%
  reshape2::melt(., id.vars=c("point", "amt_tc", "pct_tc")) %>%
  mutate(variable = case_when(variable == "contig" ~ "(c) Forest contiguity (r = 0.86)", 
                              variable == "dist_from_edge" ~ "(b) Distance to closest edge (r = 0.87)",
                              variable == "edge" ~ "(a) Area of forest edge (r = -0.90)", 
                              variable == "pct_tc" ~ "Percentage forest (lscape)", 
                              variable == "pct_core_area" ~ "(d) Percentage core area (r = 0.90)")) %>%
  ggplot(aes(pct_tc, value)) + geom_point(alpha=.5) +
  facet_wrap(~variable, ncol=2, scales="free") +
  theme_bw() +
  theme(aspect.ratio = 0.8, strip.background = element_blank(), 
        strip.text = element_text(face="bold", hjust=0), 
        panel.grid=element_blank(), 
        axis.text = element_text(colour="black")) +
  labs(x = "Percentage forest", y = "Value") 

ggsave("figures/landscape forest metrics.png", width=8.05, height=6.1, dpi=400)


## Point elevations & configuration ----
points_sp <- st_as_sf(df_ptInfo, coords=c("lon", "lat"), crs="epsg:4326") %>%
    st_transform(., crs="epsg:3117") %>%
    filter(ele_jaxa >= 880 & site_code_2 != "CC") %>%
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
    guides(fill="none") +
    labs(x="", y="")

p2 <- ggplot(points_sp, aes(rot, ele_jaxa, fill=ele_jaxa)) + 
    geom_jitter(width=.05, height=0, pch=21) + 
    scale_fill_gradientn(colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                                     "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43",
                                     "#D73027")) +
    scale_y_continuous(breaks=seq(0, 4000, 1000)) +
    # scale_x_discrete(expand=c(0,0.3)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5), 
          aspect.ratio=.8, 
          axis.text=element_text(colour="black"), 
          panel.grid.major = element_line(linetype = 1, colour="grey95")) +
    labs(y = "Elevation (m)", x="Site") +
    guides(fill="none")

ggsave("figures/points_by_rotation.png", p1, width=150, height=200, units="mm")
ggsave("figures/elevation_by_site.png", p2, width=205, height=155, units="mm")

## Pairwise distances ----
df_i <- points_sp
id_pt <- points_sp$point
id_cl <- points_sp$cluster
id_rot <- df_i$rotation
names(id_cl) <- id_pt
names(id_rot) <- id_pt

d <- dist(st_coordinates(points_sp)) %>%
    as.numeric
pt_combn <- combn(id_pt, m=2) 

## all pairwise distance combinations
dists_all <- tibble(x = pt_combn[1,], y = pt_combn[2,]) %>%
    mutate_all(as.character) %>%
    mutate(x_cl = recode(x, !!!id_cl), 
           y_cl = recode(y, !!!id_cl), 
           within_cl = x_cl == y_cl,
           x_rot = recode(x, !!!id_rot), 
           y_rot = recode(y, !!!id_rot), 
           within_rot = x_rot == y_rot,
           class = case_when(x_cl == y_cl ~ "within-cl", 
                             x_rot == y_rot ~ "within-rot"),
           class = ifelse(is.na(class), "between-rot", class),
           class = factor(class, levels=c("within-cl", "within-rot", "between-rot")),
           d = d)   
# remove between site distances
dists <- dists_all %>% filter(class != "between-rot")

# cols <- #colorRampPalette(
#     c("#4575B4", "#D73027", "#ABD9E9")#, #"#E0F3F8", "#FFFFBF", 
#                            #"#FEE090", 
#                            #"#FDAE61", 
#                            #"#F46D43", 
#       "#D73027")#)(3)
dist_lims <- c(0, max(dists$d))
plot_cl_dist_1 <- ggplot(dists, aes(d, fill=class)) + 
    geom_histogram(position="identity", alpha=.5, binwidth=50, boundary=0) +
    labs(x="Distance (m)", y="Frequency", fill="Within cluster") +
    scale_x_continuous(expand=c(0.01, 0), breaks=seq(0, 5000, 500)) +
    scale_y_continuous(expand=c(0.01, 0)) +
    scale_fill_manual(values = cols) +
    guides(fill=F) +
    theme_classic() +
    coord_cartesian(xlim=dist_lims) + 
    theme(aspect.ratio = .6, 
          # axis.title.x = element_blank(),
          # axis.text.x = element_blank(), 
          axis.line.y = element_line(size=1))

plot_cl_dist_2 <- ggplot(dists, aes(d, col=class)) + 
    geom_vline(xintercept=c(200, 500), lty=2, col="grey85") +
    stat_ecdf() +
    # geom_histogram(position="identity", alpha=.5, binwidth=100, boundary=0) +
    labs(x="Distance (m)", y="Cumulative proportion", fill="Within cluster") +
    scale_x_continuous(expand=c(0.01, 0), breaks=seq(0, 5000, 500)) +
    scale_y_continuous(expand=c(0.01, 0)) +
    scale_colour_manual(values=cols) +
    guides(col=F) +
    theme_classic() +
    coord_cartesian(xlim=dist_lims) + 
    theme(aspect.ratio = .6, 
          axis.line.y = element_line(size=1),
          # axis.title.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey95"))

plot_cl_dist_3 <- ggplot(dists_all, aes(log(d), fill=class)) + 
    geom_histogram(position="identity", alpha=.5, binwidth=.1, boundary=0) +
    labs(x="Distance (m)", y="Frequency", fill="Within cluster") +
    scale_x_continuous(expand=c(0.01, 0), #limits=log(dist_lims+c(150,200)), 
                       breaks=log(seq(0, 5000, 500)), labels=seq(0, 5000, 500)) +
    scale_y_continuous(expand=c(0.01, 0)) +
    scale_fill_manual(values = cols) +
    guides(fill=F) +
    theme_classic() +
    theme(aspect.ratio = .6, 
          # axis.title.x = element_blank(),
          # axis.text.x = element_blank(), 
          axis.line.y = element_line(size=1))

## within cluster distances ----
(mean_dists <- dists %>% group_by(within_cl) %>%
    summarise(sum(d >= 190 & d <= 500)/n(), 
              median(d)))

(mean_dists_site <- dists_all %>% group_by(within_rot) %>%
    summarise(sum(d > 500)/n(), 
              median(d)))

analysis_df %>%
  select(point, n_visit) %>%
  unique %>%
  summarise(n = sum(n_visit != 4), 
            p = n/n())

# Summary stats ----
# number of point visits 
analysis_df %>% 
  select(point, n_visit) %>%
  unique %>%
  summarise(sum(n_visit))
# number of point:species combinations
sum(analysis_df$n_visit)
# number of points
nrow(analysis_df)
# number of point:species detections
sum(analysis_df$Q)

# total number of detections
analysis_df %>% 
  summarise_at(paste0("d", 1:4), function(x) sum(x==1)) %>%
  rowSums()

## Species:point detections ----
points_by_species <- analysis_df %>% 
    filter(Q == 1) %>%
    group_by(species) %>%
    summarise(n_Q = n()) 

# 52% of species detected on fewer than 4 points
points_by_species %>% summarise(sum(n_Q < 4)/n())
points_by_species %>% summarise(sum(n_Q < 3)/n())
points_by_species %>% summarise(sum(n_Q < 2)/n())
points_by_species %>% summarise(sum(n_Q == 1)/n())

plot_points_by_species1 <- ggplot(points_by_species, aes(n_Q)) + #y=cumsum(..count..)/sum(..count..))) +
    # geom_histogram(boundary=0.5, binwidth=1, fill="grey90", colour="black") +
  stat_ecdf() +
    scale_x_continuous(breaks=c(2, 4, 8, 16, 32, 64, 128)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    theme(axis.text = element_text(colour="black"),
          axis.line.y = element_line(colour="black", size=1),
          panel.grid.major.y = element_line(colour="grey95"),
          aspect.ratio = .5) +
    labs(x = "Number of point-detections per species",
         y = "Cumulative proportion")

plot_points_by_species2 <- ggplot(points_by_species, aes(n_Q, y=..count../sum(..count..))) + 
    geom_histogram(boundary=0.5, binwidth=1, fill="grey90", colour="black") +
    scale_x_continuous(breaks=c(2, 4, 8, 16, 32, 64, 128)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0, .3, .1)) +
    theme_classic() +
    theme(axis.text = element_text(colour="black"), 
          axis.line.y = element_line(colour="black", size=1),
          panel.grid.major.y = element_line(colour="grey95"), 
          aspect.ratio = .5, 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank()) +
    labs(x = "Number of point-detections per species", 
         y = "Proportion") 
plot_points_by_species_both <- egg::ggarrange(plot_points_by_species2, 
                                              plot_points_by_species1, 
                                              ncol=1)
ggsave("figures/density_points_by_species.png",
       plot = plot_points_by_species_both, 
       height=150, width=150, units="mm", dpi = 150)

## Elevation histogram ----
eleRange_lims <- analysis_df$ele_sc %>% range

scaled_eles <- analysis_df %>% filter(Q==1) %>% pull(ele_sc)
sum(abs(scaled_eles) <=1)/length(scaled_eles)
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

props <- analysis_df %>% 
  mutate(ele_f = cut(ele_sc, seq(-50, 50, .2), 
                     labels = seq(-49.9, 49.9, .2))) %>%
    group_by(Q, ele_f) %>%
    summarise(N = n()) %>%
    group_by(ele_f) %>% 
    # summarise(N_tot = sum(N)) %>%
    summarise(N_rel = N[2]/sum(N)) %>%
    mutate(ele_f = as.numeric(as.character(ele_f)))

analysis_df %>% filter(Q==1) %>% summarise(median(ele_sc))
 
eleHist_2 <- ggplot(props, aes(ele_f, N_rel)) + 
    # geom_col(col="black", width = .1,) +
    stat_summary(fun=sum,geom="bar", colour="black", fill="grey50", width=.2) +
    scale_y_continuous(expand = expansion(0), breaks=seq(0, .2, .04)) +
    scale_x_continuous(expand = expansion(0), breaks=c(-5, -2, -1, 0, 1, 2, 5), 
                       limits = eleRange_lims) +
    theme(panel.background = element_blank(), 
          axis.line.x = element_line(), 
          axis.line.y = element_line(),
          panel.grid = element_blank(),
          axis.text=element_text(colour="black")) +
    labs(x="Range position", y="Proportion") #+
    # geom_rangeframe(data=data.frame(x=c(-5,5), y=c(0,.2)), aes(x, y))

p_ele_both <- egg::ggarrange(eleHist_1, eleHist_2, heights=c(1,.3))
ggsave("figures/elevational_histograms.png", p_ele_both, 
       width=200, height=120, units="mm", dpi=200)

## Elevation difference ----
# 1695 species:point combinations within published bounds
analysis_df %>% filter(Q==1) %>% filter(abs(ele_sc) <=1) %>% nrow
analysis_df %>% filter(Q==1) %>% filter(abs(ele_sc)  > 1) %>% nrow

1249/(1249+75)

upr_discrepancy <- analysis_df_full %>%
  mutate(upr_both = mean(c(upr_McMullan, upr_Quinones)), 
         upr_both = ifelse(is.na(upr_both), upr_Quinones, upr_both), 
         upr_diff = ele - upr_both) %>% 
  filter(ele > upr_both, Q==1) 

lwr_discrepancy <- analysis_df_full %>%
  mutate(lwr_both = mean(c(lwr_McMullan, lwr_Quinones)), 
         lwr_both = ifelse(is.na(lwr_both), lwr_Quinones, lwr_both), 
         lwr_diff = ele - lwr_both) %>% 
  filter(ele < lwr_both, Q==1) %>%
  mutate(lwr_diff = abs(lwr_diff))

lwr_discrepancy %>% filter(lwr_diff > 500)  
upr_discrepancy %>% filter(upr_diff > 500)  

plot_disc1 <- ggplot(upr_discrepancy, aes(upr_diff)) + 
  geom_histogram(boundary=0, binwidth=25, fill="grey90", col="black") +
  scale_y_continuous(breaks=seq(0, 100, 1)) +
  scale_x_continuous(limits=c(0, 800)) +
  scale_fill_manual(values=c("grey50")) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(),
        plot.title = element_text(size=12),
        panel.grid = element_blank(),
        axis.text=element_text(colour="black")) + 
  labs(y = "Frequency", 
       title = "(a) Upper-range margin")

plot_disc2 <- ggplot(lwr_discrepancy, aes(lwr_diff)) + 
  geom_histogram(boundary=0, binwidth=25, fill="grey90", col="black") +
  scale_y_continuous(breaks=seq(0, 100, 1)) +
  scale_x_continuous(limits=c(0, 800)) +
  scale_fill_manual(values=c("grey50")) +
  theme_bw() +
  theme(panel.background = element_blank(), 
        axis.line.y = element_line(),
        panel.grid = element_blank(),
        plot.title = element_text(size=12),
        axis.text=element_text(colour="black")) + 
  guides(fill=FALSE) +
  labs(x = "Distance from range margin (m)", 
       y = "Frequency", 
       title = "(b) Lower-range margin")

pboth_disc <- egg::ggarrange(plot_disc1, plot_disc2)
png("figures/elevational_limit_discrepancy.png", width=150, height=150, units="mm", res=100)
pboth_disc
dev.off()
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
png("figures/rangeMargin_comparison.png", height=100, width=150, units="mm", res=100) 
rangeMar_both
dev.off()


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


## analysis 
pt_visit <- analysis_df %>% 
  select(n_visit, point) %>% 
  unique
pt_visit %>% filter(n_visit == 2)
pt_visit %>% filter(n_visit == 3)
