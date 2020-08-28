## Create forest species subset 
library(dplyr); library(ggplot2)
parker_lookup <- readRDS("data/Parker_Stotz_Fitzpatrick_1996/parker_lookup_subset.rds")
adata <- read.csv("data/Parker_Stotz_Fitzpatrick_1996/databases/adata.csv") %>%
    as_tibble %>%
    mutate(parker = paste(GENUS, SPECIES))

species_lookup <- read.csv("data/initial_species_list.csv") %>%
    left_join(., parker_lookup) %>%
    as_tibble %>%
    mutate(Species_HBW = gsub(" ", "_", HBW)) %>%
    select(Species_HBW, eBird, HBW, Donegan, Pulido, eltontraits, parker)

df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds") 

# Drops 1 species: Columba livia
df_species <- df_bird %>% 
    group_by(Species_HBW) %>%
    summarise(n_tot = n(), n_pt = length(unique(Point))) %>%
    select(Species_HBW, n_tot, n_pt) %>%
    arrange(desc(n_pt)) %>%
    left_join(., species_lookup) %>%
    left_join(., adata) 

dropped_species <- df_species %>%
    filter(is.na(parker))
dropped_species

df_summ <- df_bird %>% 
    group_by(Species_HBW) %>%
    mutate(n_tot = n()) %>%
    group_by(Species_HBW, n_tot, Point) %>%
    summarise(Q = 1) %>%
    group_by(Species_HBW, n_tot) %>%
    summarise(n_pt = sum(Q)) %>%
    arrange(desc(n_pt))

library(ggplot2)
ggplot(df_summ %>% filter(n_pt<10), aes(n_pt, n_tot)) + geom_jitter(height=.2, width=.2) + 
    scale_y_continuous(breaks=seq(0, 20, 2)) +
    scale_x_continuous(breaks = seq(0, 20, 2)) + 
    coord_equal()


df_habcat <- df_species %>% select(Species_HBW, n_pt, n_tot, HAB1:HAB7) %>%
    reshape2::melt(., id.vars=c("Species_HBW", "n_pt", "n_tot")) %>%
    select(Species_HBW, n_pt, n_tot, order=variable, habitat=value)

df_habcat_summ <- df_habcat %>% group_by(Species_HBW, n_pt, n_tot) %>%
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
