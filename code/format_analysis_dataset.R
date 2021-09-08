# Formatting dataset for running analysis and passing to stan
# Notes:
# (1) analysis_df is the full final analysis dataframe
# (2) stan_data is the same data, repackaged for passing to stan
# (3) gets called from model run scripts

# packages
library("dplyr"); library("reshape2"); library("cmdstanr"); library("posterior")

# functions to scale and unscale elevational range
# rewrite these to be in terms of Upper and Lower
range_scaling <- function(ele, mid, scale) {
    (ele - mid)/(scale)
}
range_unscaling <- function(ele_scaled, mid, scale) {
    ele_scaled *scale + mid
}

## Read in datasets ----
#* species lookup table ----
species_lookup <- read.csv("data/initial_species_list.csv", as.is=T) 
colnames(species_lookup) <- tolower(colnames(species_lookup))
species_lookup <- species_lookup %>%
    mutate(species_clements = gsub(" ", "_", ebird), 
           species_hbw = gsub(" ", "_", hbw)) %>%
    select(-x)

#* birdlife habitat specialism ----
# cross-reference with forest specialism subset
species_forest <- readRDS("data/birdlife_traits_df.rds") %>%
    mutate(species = gsub(" ", "_", species))

#* detection dataset ----
df_bird <- readRDS("data/analysisDataset_EasternCordillera.rds") %>%
    dplyr::select(-site) %>%
    left_join(., species_lookup, by="species_hbw") %>%
    mutate(species_clements = species_clements.y) %>%
    dplyr::select(-species_clements.x, -species_clements.y)

#* elevational ranges (Quinones) ----
ele_Quinones <- read.csv("data/elevational_ranges_Quinones.csv", as.is=T) %>%
    mutate(species_clements = paste0(genus, "_", species)) %>%
    mutate(lower = ifelse(is.na(lower), 0, lower)) %>%
    left_join(species_lookup, .) %>%
    mutate(species = species_hbw) %>% 
    as_tibble %>%
    rename(lwr_Quinones = lower, 
           upr_Quinones = upper)

#* Elevational ranges (Mcmullan) ----
ele_McMullan <- read.csv("data/Jacob_dropbox/Bird_elevations_initial.csv") %>%
    as.data.frame %>%
    mutate(species = gsub(" ", "_", Scientific)) %>%
    as_tibble %>%
    select(species, lwr_McMullan=eMin, upr_McMullan=eMax) 

#* forest cover ----
# 50% forest cover threshold, and remove loss 
forest_cover <- readRDS("data/hansen_500m_buffer.rds") %>%
    mutate(tc_2000 = treecover2000 >= 50, 
           tc_2018 = as.numeric(tc_2000 - (lossyear != 0) == 1)) %>%
    group_by(point) %>%
    summarise(amt_tc = sum(tc_2018), pct_tc = amt_tc/n())

#* Point info ----
df_ptInfo <- readRDS("data/point_ele_EasternCordillera.rds") %>%
    select(-geometry) %>%
    as.data.frame() %>%
    mutate(ele = ele_jaxa) %>%
    left_join(., forest_cover, by="point") %>% 
    filter(forest == 1) %>%
    mutate(habitat = pct_tc, 
           rotation = site_code_2)

#* migration timings ----
migrant_dates <- read.csv("data/Jacob_dropbox/migratory.csv") %>%
    as_tibble %>%
    select(species = 1, status:note) %>%
    mutate(species = gsub(" ", "_", species)) %>%
    filter(start1 != "") %>%
    mutate(start1 = lubridate::yday(as.Date(start1, format="%b_%d")), 
           end1 = lubridate::yday(as.Date(end1, format="%b_%d")))

#* Calculate average ele limits ----
ele_both <- full_join(ele_Quinones, ele_McMullan) %>%
    mutate(midpoint_Quinones = (upr_Quinones + lwr_Quinones)/2, 
           scale_Quinones = (upr_Quinones - lwr_Quinones)/2, 
           midpoint_McMullan = (upr_McMullan + lwr_McMullan)/2, 
           scale_McMullan = (upr_McMullan - lwr_McMullan)/2, 
           midpoint_both = case_when(is.na(midpoint_McMullan) ~  midpoint_Quinones,
                                     is.na(midpoint_Quinones) ~ midpoint_McMullan, 
                                     T ~ (midpoint_McMullan + midpoint_Quinones)/2),
           scale_both = case_when(is.na(scale_McMullan) ~ scale_Quinones,
                                  is.na(scale_Quinones) ~ scale_McMullan, 
                                  T ~ (scale_McMullan + scale_Quinones)/2), 
           upr_both = midpoint_both + scale_both, 
           lwr_both = midpoint_both - scale_both)


## Get full df and zf version ----
# filtering for forest points (note this df still retains 0-bird visits)
# note: !duplicated() is to drop 21 point-visits where a species is recorded 
# twice. Also remove 3 migrant species.
forest_specialist_subset <- df_bird %>% 
    left_join(., df_ptInfo) %>% 
    group_by(species_hbw, point) %>% 
    slice(1) %>%
    group_by(species_hbw, forest) %>%
    summarise(n = n()) %>%
    group_by(species_hbw) %>% mutate(n_tot = sum(n)) %>%
    filter(forest == 1) %>%
    mutate(prop = n/n_tot, species=species_hbw) %>%
    left_join(., species_forest %>% rename(species_hbw=species))

df_det <- df_bird %>% 
    left_join(., df_ptInfo) %>%
    filter(!(species_hbw %in% migrant_dates$species)) %>%
    filter(forest == 1) %>%
    filter(ele >= 880) %>%
    group_by(point, visit) %>%
    filter(!duplicated(species_hbw)) %>%
    ungroup %>%
    mutate(det = 1) %>%
    dplyr::select(site, rotation, cluster, 
                  species = species_hbw, point, visit, date, 
                  time, det, observer) 

# get all point:visit combinations
point_visits <- df_det %>% select(point, visit) %>% unique

#* subset for forest specialists ----
df_subset <- df_det %>%
    left_join(., forest_specialist_subset) %>% 
    left_join(., ele_both) %>%
    left_join(., species_forest) %>%
    group_by(species) %>%
    mutate(sumQ = length(unique(point))) %>%
    ungroup %>%
    filter(forest_dep %in% c("High"))

unique_visits <- point_visits %>%
    left_join(., df_bird[c("point", "visit", "date", "time", "observer", "no_birds")] %>% 
                  unique) %>%
    left_join(., df_ptInfo)

unique_species <- unique(df_subset$species)
n_species <- length(unique_species)

#* zero-fill dataframe ----
zf_subset <- replicate(n_species, unique_visits, simplify=F) %>%
    bind_rows(., .id="species") %>%
    mutate(species = unique_species[as.numeric(species)]) %>%
    # retain rows that aren't in df_det, specify detection is 0, then add det
    anti_join(., df_subset) %>%
    mutate(det = 0) %>%
    bind_rows(., df_subset) %>%
    # remove the 'NA' species (created by points that were 0-birds)
    filter(!is.na(species)) %>%
    mutate(time_sc = (time-median(time))/sd(time), 
           obsvr_num = as.numeric(as.factor(observer)), 
           site_obs = paste0(rotation, observer), 
           site_obs_num = as.integer(as.factor(site_obs))) %>%
    group_by(point) %>%
    mutate(visit_new = as.integer(as.factor(visit))) %>%
    ungroup 

observer_ptLevel <- zf_subset %>%
    group_by(point) %>%
    summarise(obsvr_point = names(sort(table(observer), TRUE)[1]))

# check: zf subset has the combinations that it should..
# zf_subset %>%
#     group_by(point) %>%
#     summarise(length(unique(observer)), names(sort(table(observer), TRUE)[1])) %>%
#     View

# check: number of species & point:visit combinations is the number of rows
# length(unique(zf_subset$species)) *
#     nrow(unique(zf_subset[,c("point", "visit")])) == nrow(zf_subset)

# Format for modelling ----
det <- dcast(zf_subset, point + species ~ visit_new, value.var="det") %>%
    rename(d1 = `1`, d2 = `2`, d3=`3`, d4 = `4`)
vis <- dcast(zf_subset, point + species ~ visit_new, value.var="time_sc") %>%
    rename(v1 = `1`, v2 = `2`, v3=`3`, v4 = `4`)
obsvr <- dcast(zf_subset, point + species ~ visit_new, value.var="observer") %>%
    rename(o1 = `1`, o2 = `2`, o3=`3`, o4 = `4`)

id_site_obs <- dcast(zf_subset, point + species ~ visit_new, value.var="site_obs_num") %>%
    rename(site_obs1 = `1`, site_obs2 = `2`, site_obs3=`3`, site_obs4 = `4`)

id_obsvr_JS <- obsvr %>% 
    mutate_at(c("o1", "o2", "o3", "o4"), function(x) as.integer(x == "JS")) %>%
    rename(JS1 = o1, JS2 = o2, JS3 = o3, JS4 = o4)

id_obsvr_DE <- obsvr %>% 
    mutate_at(c("o1", "o2", "o3", "o4"), function(x) as.integer(x == "DPE")) %>%
    rename(DE1 = o1, DE2 = o2, DE3 = o3, DE4 = o4)

# Set NA to be -99
det[is.na(det)] <- -99
vis[is.na(vis)] <- -99
obsvr[is.na(obsvr)] <- -99
id_obsvr_JS[is.na(id_obsvr_JS)] <- -99
id_obsvr_DE[is.na(id_obsvr_DE)] <- -99
id_site_obs[is.na(id_site_obs)] <- 1 # needs to be 1 b/c of model coding

det_vis <- left_join(det, vis) %>%
    left_join(., id_obsvr_JS) %>%
    left_join(., id_obsvr_DE) 

analysis_df <- left_join(det_vis, df_ptInfo) %>%
    left_join(., ele_both) %>%
    left_join(., obsvr) %>%
    left_join(., observer_ptLevel) %>%
    left_join(., species_forest[c("species", "forest_dep")]) %>%
    ungroup %>%
    filter(!species %in% c("Vireo_chivi", "Vireolanius_eximius")) %>%
    mutate(ele = round(ele, 0),
           ele_std = (ele-mean(ele))/sd(ele),
           ele_sc = range_scaling(ele, midpoint_both, scale_both), 
           habitat_sc = (habitat - mean(habitat))/sd(habitat),
           ele_median_sc = (midpoint_both - mean(midpoint_both))/sd(midpoint_both),
           ele_midpoint_cut = cut(midpoint_both, seq(200, 4000, 200)), 
           ele_midpoint_f = as.integer(ele_midpoint_cut), 
           ele_band = cut(ele, seq(min(ele), max(ele), len=5), include.lowest = T),
           distance_lwr = ifelse(ele <= midpoint_both, ele - lwr_both, upr_both - ele), 
           distance_upr = ifelse(ele > midpoint_both, 0, upr_both - ele),
           n_visit = rowSums((select(., d1:d4) != -99)), 
           Q = as.integer(apply(select(., d1:d4), 1, function(x) any(x == 1)))) %>%
    mutate(id_sp = as.numeric(as.factor(species)),
           id_pt = as.numeric(as.factor(point)),
           id_site_sp = paste0(rotation, species),
           id_site_sp_num = as.numeric(as.factor(id_site_sp)), 
           id_cl_sp = paste0(cluster, species),
           id_cl_sp_num = as.numeric(as.factor(id_cl_sp)), 
           id_obs_sp = paste0(obsvr_point, species), 
           id_obs_sp_num = as.numeric(as.factor(id_obs_sp)), 
           id_ele_band = as.integer(ele_band)) %>%
    select(-donegan, -donegan2, -ebird2, -hbw, -ebird, -pulido, -eltontraits, 
           -english, -genus, -is_na, -ele_gps, -ele_srtm) %>%
    select(rotation, cluster, point, species, d1:d4, Q, everything()) 

print(paste0("nrow: ", nrow(analysis_df)))
print(paste0("nspecies: ", length(unique(analysis_df$species))))

# create lookup relating species number to name
lookup <- tibble(species = as.factor(det$species), 
                 id_sp = as.numeric(as.factor(det$species))) %>% unique %>%
    mutate(id_sp = as.numeric(id_sp))

stan_data <- list(n_visit = analysis_df$n_visit,
                  n_species = length(unique(analysis_df$species)),
                  n_points = length(unique(analysis_df$point)),
                  n_tot = nrow(analysis_df),
                  n_ele_midpoint_f = length(unique(analysis_df$ele_midpoint_f)),
                  id_sp = analysis_df$id_sp,
                  id_pt = analysis_df$id_pt,
                  id_cl_sp = analysis_df$id_cl_sp_num,
                  id_site_sp = analysis_df$id_site_sp_num,
                  id_obs_sp = analysis_df$id_obs_sp_num,
                  id_ele_band = analysis_df$id_ele_band,
                  id_obs_JS = analysis_df %>% select(JS1:JS4), 
                  id_obs_DE = analysis_df %>% select(DE1:DE4), 
                  n_cluster_species = max(analysis_df$id_cl_sp_num),
                  n_site_species = max(analysis_df$id_site_sp_num),
                  n_obs_species = max(analysis_df$id_obs_sp_num),
                  Q = analysis_df$Q,
                  ele_midpoint = analysis_df$ele_median_sc,
                  ele_midpoint_f = analysis_df$ele_midpoint_f,
                  ele = analysis_df$ele_std,
                  distance_lwr = analysis_df$distance_lwr,
                  distance_upr = analysis_df$distance_upr,
                  time = analysis_df %>% select(v1:v4),
                  det_data = analysis_df %>% select(d1:d4),
                  habitat = analysis_df$habitat_sc,
                  range_pos = analysis_df$ele_sc,
                  grainsize=1)

# Save data ----
saveRDS(stan_data, "output/stan_data.rds")
saveRDS(analysis_df, "output/analysis_df_birdlife.rds") 
