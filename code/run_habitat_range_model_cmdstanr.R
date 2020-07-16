## run elevation model

## Housekeeping ----
library("dplyr"); library("reshape2"); library("cmdstanr"); library("posterior")

# functions to scale and unscale elevational range
range_scaling <- function(ele, mid, scale) {
    (ele - mid)/(scale)
}
range_unscaling <- function(ele_scaled, mid, scale) {
    ele_scaled *scale + mid
}


## Read in datasets ----
species_lookup <- read.csv("data/initial_species_list.csv", as.is=T) %>%
    mutate(Species_Clements = gsub(" ", "_", eBird), 
           Species_HBW = gsub(" ", "_", HBW)) %>%
    select(-X)

df_bird <- readRDS("data/analysisDataset_EasternCordillera.rds") %>%
    select(-Site) %>%
    left_join(., species_lookup, by="Species_HBW") %>%
    mutate(Species_Clements = Species_Clements.y) %>%
    select(-Species_Clements.x, -Species_Clements.y)

#* forest and landscape metrics ----
forest_cover <- readRDS("data/hansen_500m_buffer.rds") %>%
    mutate(tc_2000 = treecover2000 > 50, 
           tc_2018 = as.numeric(tc_2000 - (lossyear != 0) == 1)) %>%
    group_by(point_id) %>%
    summarise(amt_tc = sum(tc_2018), pct_tc = amt_tc/n())

#* elevational ranges (Quinones) ----
ele_range <- read.csv("data/elevational_ranges_Quinones.csv", as.is=T) %>%
    mutate(Species_Clements = paste0(Genus, "_", Species)) %>%
    mutate(Lower = ifelse(is.na(Lower), 0, Lower)) %>%
    left_join(species_lookup, .) %>%
    mutate(Species = Species_HBW) %>% 
    as_tibble

#* landscape vars----
landscape_vars <- readRDS("data_landscapeMetrics_draft_v1.rds") %>%
    select(Point, Habitat) %>%
    unique

#* point info ----
df_ptInfo <- read.csv("data/CO_sampling_points_metafile_ver2.csv") %>%
    as_tibble() %>%
    left_join(., forest_cover) %>%
    rename(Point = point_id) %>%
    left_join(., landscape_vars) %>%
    mutate(site_cluster = paste0(site, "_",cluster)) %>%
    dplyr::select(Point, lat, long, ele=ALOSelev, amt_tc, pct_tc, Site=site_cluster, 
                  Habitat)

# ggplot(df_ptInfo, aes(ele, pct_tc, col=Habitat)) + geom_point()
# df_ptInfo %>% 
#     filter(Habitat == "Forest") %>%
#     summarise(quantile(pct_tc, c(.1,.9), na.rm=T))
# (pi*100^2)/(pi * 500^2)

## Get full df and zf version ----
# filtering for forest points (note this df still retains 0-bird visits)
df_det <- df_bird %>% 
    left_join(., df_ptInfo) %>%
    filter(Habitat == "Forest") %>%
    group_by(Point, Visit) %>%
    filter(!duplicated(Species_HBW)) %>%
    ungroup %>%
    mutate(det = 1) %>%
    dplyr::select(Site, Species = Species_HBW, Point, Visit, Time, det, Observer) 

# get all point:visit combinations
point_visits <- df_det %>% select(Point, Visit) %>% unique

# note: 72 point:visit combinations that have duplicated species-- doesn't matter
# for this modelling exercise, but tidy up regardless
# det %>% group_by(Point, Visit) %>%
#     filter(duplicated(Species))

uniqueVisits <- point_visits %>%
    left_join(., df_ptInfo) %>%
    left_join(., df_bird[c("Point", "Visit", "Time", "Observer", "no_birds")] %>% unique) %>%
    mutate(Time = ifelse(is.na(Time), median(Time, na.rm=T), Time))

uniqueSpecies <- unique(df_det$Species)
n_species <- length(uniqueSpecies)

# zero-filled dataframe
zf_det <- replicate(n_species, uniqueVisits, simplify=F) %>%
    bind_rows(., .id="Species") %>%
    mutate(Species = uniqueSpecies[as.numeric(Species)]) %>%
    # retain rows that aren't in df_det, specify detection is 0, then add det
    anti_join(., df_det) %>%
    mutate(det = 0) %>%
    bind_rows(., df_det) %>%
    # remove Species==NA (from points that were 0-birds)
    filter(!is.na(Species))

#* Take subset ----
# Filter for species that are (a) montane, (b) forest specialisits, and (c) 
# present on at least 3 points
species_forest <- readRDS("data/forest_subset.rds") %>%
    left_join(., ele_range) %>%
    filter(n_pt >= 15) %>% 
    mutate(Species = Species_HBW)

zf_subset <- filter(zf_det, Species %in% species_forest$Species) %>%
    mutate(Time_sc = (Time-median(Time))/sd(Time), 
           obsvr_num = as.numeric(as.factor(Observer))) %>%
    group_by(Point) %>%
    mutate(Visit_new = as.integer(as.factor(Visit))) %>%
    ungroup

observer_lookup <- zf_subset %>% select(Observer, obsvr_num) %>% unique

# spread 
det <- dcast(zf_subset, Point + Species ~ Visit_new, value.var="det")
vis <- dcast(zf_subset, Point + Species ~ Visit_new, value.var="Time_sc")
obsvr <- dcast(zf_subset, Point + Species ~ Visit_new, value.var="obsvr_num")
# not sure if Stan actually requires this, but have set NA to be -99 anyway
det[is.na(det)] <- -99
vis[is.na(vis)] <- -99
obsvr[is.na(obsvr)] <- -99

eles <- left_join(det, df_ptInfo) %>%
    left_join(., ele_range) %>%
    mutate(midpoint = (Lower + Upper)/2, 
           scale = (Upper - Lower)/2, 
           ele_sc = range_scaling(ele, midpoint, scale), 
           habitat_sc = (pct_tc - mean(pct_tc))/sd(pct_tc),
           habitat_ele = ele_sc * habitat_sc,
           id_site_sp = interaction(substr(Point, 1, 2), Species), 
           id_site_sp_num = as.numeric(id_site_sp), 
           Q = as.numeric(rowSums(select(., `1`:`4`), na.rm = T)>0), 
           n_visit = rowSums(!is.na(select(., `1`:`4`)))) %>%
    # tidy up a bit (drop redundant cols) 
    select(-Donegan, -Donegan2, -eBird2, -HBW, -eBird, -Pulido, -eltontraits, 
           -English..part.filled., -Genus, -is_NA, -extra, -Alternative, 
           -Has_alternative)

# check there aren't NAs in dataframe
missing_ele <- eles %>% filter(is.na(Upper) | is.na(Lower)) %>% pull(Species) %>% unique
if(length(missing_ele) > 0) stop("some species are missing elevation")

lookup <- tibble(Species = as.factor(det$Species), 
                 id_sp = as.numeric(as.factor(det$Species))) %>% unique %>%
    mutate(id_sp = as.numeric(id_sp))

stan_data <- list(n_visit = eles$n_visit,
                  n_species = length(unique(det$Species)),
                  n_points = length(unique(det$Point)),
                  n_tot = nrow(det),
                  id_sp = as.numeric(as.factor(det$Species)),
                  id_pt = as.numeric(as.factor(det$Point)),
                  id_site_sp = eles$id_site_sp_num,
                  id_obsvr = obsvr[,-c(1,2)],
                  n_sites_species = max(eles$id_site_sp_num),
                  Q = eles$Q,
                  vis_cov = vis[,-c(1,2)],
                  det_data = det[,-c(1,2)],
                  ele = eles$ele_sc,
                  habitat = eles$habitat_sc,
                  midpoint = eles$midpoint,
                  grainsize=1, 
                  species_lookup = lookup) # not needed for stan, but useful later 

## Run mod ----
mod <- cmdstan_model("code/stan_files/elevation_detection_habitat_asymmetric_reparam.stan", 
                     cpp_options = list(stan_threads = T))

samps <- mod$sample(data = stan_data, 
                    chains = 2, 
                    parallel_chains = 2, #n_chains, 
                    threads_per_chain = 7, #n_threads, 
                    iter_warmup = 1000, 
                    iter_sampling = 2000)

draws <- samps$draws()
# draws without junk
draws_clean_1 <- draws[,,!grepl("raw", variables(draws))]

save_time <- Sys.time()
# save files
saveRDS(samps, paste0("output/model_habitatModel_asymm_", save_time, ".rds"))
saveRDS(draws_clean_1, paste0("output/draws_habitatModel_asymm_", save_time, ".rds"))
saveRDS(stan_data, paste0("output/stanData_habitatModel_asymm_", save_time, ".rds"))
