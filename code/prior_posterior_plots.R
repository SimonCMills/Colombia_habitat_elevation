# Prior-posterior plots
fixefs <- fixef(fit, summary=F) %>%
    reshape2::melt(., id.vars="draw") %>%
    as_tibble

x <- seq(-10, 5, len=200)
coef_names <- unique(fixefs$variable)

prior_dens <- expand.grid(variable = coef_names, x = x) %>%
    mutate(y = case_when(variable == "Intercept" ~ dnorm(x, -2.5, 1.5), 
                         variable == "occ_Intercept" ~ dnorm(x, -2.5, 1.5), 
                         variable == "occ_range_pos" ~ dnorm(x, 0, 2), 
                         variable == "occ_range_pos2" ~ dnorm(x, -4, 2), 
                         
                         variable == "occ_ele_midpoint_std" ~ dnorm(x, 0, 2), 
                         variable == "occ_ele_midpoint_std:habitat_std" ~ dnorm(x, 0, 2), 
                         variable == "occ_range_pos:habitat_std" ~ dnorm(x, 0, 2), 
                         variable == "occ_habitat_std:range_pos_upr" ~ dnorm(x, 0, 2), 
                         T ~  dnorm(x, 0, .5))) %>%
    filter(y > 0.01)


ggplot(fixefs, aes(value)) + geom_density() +
  geom_ribbon(data=prior_dens, aes(x=x, ymin=0, ymax=y), alpha=.4) +
  facet_wrap(~variable, scales="free")

# sd components
vc <- VarCorr(fit, summary = F)
vc$id_obs_sp

vc2 <- lapply(vc, as_tibble)
vc_list <- list()
for(i in 1:length(vc)) {
  vc_i <- as_tibble(vc[[i]]$sd)
  names(vc_i) <- paste0(names(vc_i), "_", names(vc)[i])
  vc_list[[i]] <- vc_i
  }

sd_df <- bind_cols(vc_list)
sd_lf <- reshape2::melt(sd_df)

x <- seq(0, 5, len=200)
prior_dens <- expand.grid(variable = names(sd_df), x = x) %>%
  mutate(y = case_when(variable == "occ_Intercept_id_cl_sp" ~ dnorm(x, 0, 3)*2, 
                       variable == "occ_Intercept_id_site_sp" ~ dnorm(x, 0, 3)*2, 
                       T ~ dnorm(x, 0, 1)*2)) %>%
  filter(y > 0.01)

ggplot(sd_lf, aes(value)) + geom_density() +
  geom_ribbon(data=prior_dens, aes(x=x, ymin=0, ymax=y), alpha=.4) +
  facet_wrap(~variable, scales="free")
