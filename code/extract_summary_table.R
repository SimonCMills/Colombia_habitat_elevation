# extract summary table for SOM
library(dplyr); library(kableExtra); library(brms); library(kableExtra)

summ <- summary(fit)

RE_cleaned <- bind_rows(summ$random, .id = "id") %>%
    mutate(Name = row.names(.)) %>%
    as_tibble %>%
    select(Name, id, Mean = Estimate, Lower_CI = `l-95% CI`, Upper_CI = `u-95% CI`) %>%
    # renaming
    mutate(Class  = ifelse(grepl(".*\\(occ_", Name), 
                           "Occupancy random effect", 
                           "Detection random effect"), 
           id = gsub("^id_", "", id),
           id = case_when(id == "rotation" ~ "site", 
                          id == "species_eltontraits" ~ "phylo", 
                          TRUE ~ id),
           Name = tolower(Name), 
           Name = gsub("sd\\(", "sd_", Name),
           Name = gsub("cor\\(", "cor_", Name), 
           Name = gsub("\\).*", "", Name), 
           Name = gsub("_std|occ_", "", Name), 
           Group = ifelse(grepl("^cor_", Name), gsub("^cor_(.*)", "\\1", Name), ""),
           Group = gsub("\\,", "\\, ", Group),
           Name = ifelse(grepl("^cor_", Name), "cor", Name), 
           Name = paste0(Name, "_", id),
           ) %>%
    select(Class, Name, Group, Mean, Lower_CI, Upper_CI)

fEff_clean <- fixef(fit) %>%
    as.data.frame %>%
    mutate(Name = row.names(.)) %>%
    as_tibble %>%
    select(Name, Mean = Estimate, Lower_CI = Q2.5, Upper_CI = Q97.5) %>%
    mutate(Class = ifelse(grepl("occ_", Name), 
                          "Occupancy fixed effect", 
                          "Detection fixed effect"), 
           Name = gsub("occ_|_point|_std", "", Name), 
           Group = "")

summ_table <- bind_rows(fEff_clean, RE_cleaned) %>%
    select(Class, Name, Group, everything()) %>%
    mutate(Class = factor(Class, levels = c(sort(unique(Class))[c(3, 4, 1, 2)], ""))) %>%
    arrange(Class)

summ_kable <- summ_table %>%
    mutate_at(vars(c("Mean", "Lower_CI", "Upper_CI")), function(x) round(x, 3)) %>%
    mutate_at(vars(c("Mean", "Lower_CI", "Upper_CI")), 
              function(x) cell_spec(x, bold = ifelse(x > 0, "T", "F"))) %>%
    rename(`Lower CI` = Lower_CI, `Upper CI` = Upper_CI) %>%
    kable(escape = F, booktabs = T) %>%
    collapse_rows(columns = 1, valign = "top") %>%
    kable_styling()

save_kable(summ_kable, file = "output/model_summary_table.html")
