## Create forest species subset from the PSF database. Forest specialists are 
## defined as those species that are present in forest, but never present in 
## forest edge, secondary forest, and pasture. 

library(dplyr); library(ggplot2)

parker_lookup <- readRDS("data/Parker_Stotz_Fitzpatrick_1996/parker_lookup_subset.rds")
adata <- read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv") %>%
    as_tibble %>%
    mutate(parker = paste(GENUS, SPECIES))

cdata <- read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/cdata.csv") %>%
    as_tibble %>%
    mutate(parker = paste(GENUS, SPECIES))

parker_data <- bind_rows(adata, cdata)
dim(adata)
dim(parker_data)

species_lookup <- read.csv("data/initial_species_list.csv") %>%
    left_join(., parker_lookup) %>%
    as_tibble %>%
    mutate(species_hbw = gsub(" ", "_", HBW)) %>%
    select(species_hbw, eBird, HBW, Donegan, Pulido, eltontraits, parker)

df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds") 

# Drops 1 species: Columba livia
df_species <- df_bird %>% 
    group_by(species_hbw) %>%
    summarise(n_tot = n(), n_pt = length(unique(point))) %>%
    select(species_hbw, n_tot, n_pt) %>%
    arrange(desc(n_pt)) %>%
    left_join(., species_lookup) %>%
    left_join(., parker_data) 

dropped_species <- df_species %>%
    filter(is.na(parker))
dropped_species

df_summ <- df_bird %>% 
    group_by(species_hbw) %>%
    mutate(n_tot = n()) %>%
    group_by(species_hbw, n_tot, point) %>%
    summarise(Q = 1) %>%
    group_by(species_hbw, n_tot) %>%
    summarise(n_pt = sum(Q)) %>%
    arrange(desc(n_pt))

ggplot(df_summ %>% filter(n_pt<=20), aes(n_pt, n_tot/4/n_pt)) + geom_jitter(height=0, width=.2) + 
    scale_y_continuous("N observations/N points", breaks=seq(0, 1, .25)) +
    scale_x_continuous("N points observed", breaks = seq(0, 20, 2)) + 
    coord_fixed(5)


df_habcat <- df_species %>% select(species_hbw, n_pt, n_tot, HAB1:HAB7) %>%
    reshape2::melt(., id.vars=c("species_hbw", "n_pt", "n_tot")) %>%
    select(species_hbw, n_pt, n_tot, order=variable, habitat=value)

df_habcat_summ <- df_habcat %>% group_by(species_hbw, n_pt, n_tot) %>%
    summarise(any_forest = any(habitat %in% paste0("F", 1:14)), 
              sec_and_forest_edge = any(habitat %in% paste0("F", 15)|habitat %in% paste0("F", 1:15, "E")), 
              pasture = any(habitat %in% paste0("N", 1:15)|habitat %in% paste0("N", 1:15, "E")),
              preference = habitat[order == "HAB1"] %in% paste0("F", 1:14))

summary_table <- df_habcat_summ %>% group_by(any_forest, sec_and_forest_edge, pasture) %>% 
    summarise(n_sp = n())

write.csv(df_habcat_summ, "data/habitatprefs.csv")

species_forest <- df_habcat_summ %>%
    filter(any_forest, !sec_and_forest_edge, !pasture) %>%
    unique() %>%
    arrange(desc(n_pt))

saveRDS(species_forest, "data/forest_subset.rds")
