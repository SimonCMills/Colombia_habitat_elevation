# Prior-posterior plots

# read draws
posterior_dist <- read_cmdstan_csv("stan_out_local/model_distance_sumQ_gt1000_chains_1-202106110055-1.csv")
draws_dist <- posterior_dist$post_warmup_draws %>% subset_draws(iteration = 1:500)

posterior_ele <- read_cmdstan_csv("stan_out_local/model_ele_linear-202107060909-1-09e59c.csv")
                                   
posterior_dist <- read_cmdstan_csv("stan_out_local/model_distance_sumQ_gt1000_chains_1-202106110055-1.csv")
draws_dist <- posterior_dist$post_warmup_draws

posterior_ele <- read_cmdstan_csv("stan_out_local/model_ele_linear-202106102157-1-5f039d.csv")
draws_ele <- posterior_ele2$post_warmup_draws

# x ranges and prior densities
x_sigma <- seq(0, 5, len=100)
sigma_terms <- dnorm(x_sigma, 0, 2)*2
plot(sigma_terms ~ x_sigma)

x_pars <- seq(-5, 5, len=100)
intercept_terms <- dnorm(x_pars, 0, 1.5)
slope_terms <- dnorm(x_pars, 0, 3)
slope_ele <- dnorm(x_pars, -3, 4)
mid_terms <- dnorm(x_pars, 0, 1)

# levels (for ordering in plot)
lvls <- c("mu_mu", paste0("mu_b", 0:2), "b3", "b4", paste0("mu_b", 5:7), 
          "sigma_mu", paste0("sigma_b", 0:2), paste0("sigma_b", 5:7), 
          paste0("mu_d", 0:1), paste0("sigma_d", 0:1), 
          "sigma_cl_sp", "sigma_site_sp", "sigma_obs_sp")

# create prior dataframe
priors <- expand.grid(type = c("intercept", "slope", "slope_ele", "mid0", "sigma"), id = 1:100) %>%
    mutate(x = case_when(type %in% c("intercept", "slope", "mid0") ~ x_pars[id], 
                         type == "slope_ele" ~ seq(-7.5, 5, len=100)[id],
                         type %in% c("sigma") ~ x_sigma[id]), 
           y = case_when(type %in% c("intercept") ~ intercept_terms[id],
                         type %in% c("slope") ~ slope_terms[id], 
                         type %in% c("slope_ele") ~ slope_ele[id], 
                         type %in% c("mid0") ~ mid_terms[id],
                         type %in% c("sigma") ~ sigma_terms[id]))

# create posterior data frame
posterior <- draws_ele %>%
    subset_draws("mu|sigma|b4|b5", regex=T) %>%
    as_draws_df() %>% 
    reshape2::melt() %>%
    as_tibble %>%
    filter(!(variable %in% c(".chain", ".iteration", ".draw"))) %>%
    mutate(type = case_when(variable %in% c("mu_b0", "mu_d0") ~ "intercept", 
                            variable %in% c("mu_b1", "mu_b3", "mu_b3_upr", "b4", "b5", "mu_d1") ~ "slope",
                            variable %in% c("mu_b2") ~ "slope_ele",
                            variable %in% c("mu_mid0") ~ "mid0", 
                            grepl("sigma", variable) ~ "sigma"), 
           variable2 = case_when(variable == "b4" ~ "b3", 
                                 variable == "b5" ~ "b4", 
                                 variable == "mu_b3" ~ "mu_b7",
                                 variable == "sigma_b3" ~ "sigma_b7",
                                 variable == "mu_mid0" ~ "mu_mu",
                                 variable == "sigma_mid0" ~ "sigma_mu",
                                 TRUE ~ as.character(variable)), 
           variable2 = factor(variable2, levels=lvls))

posterior2 <- draws_dist %>%
    subset_draws("mu|sigma|b4|b5", regex=T) %>%
    as_draws_df() %>% 
    reshape2::melt() %>%
    as_tibble %>%
    filter(!(variable %in% c(".chain", ".iteration", ".draw"))) %>%
    mutate(type = case_when(variable %in% c("mu_b0", "mu_d0") ~ "intercept", 
                            variable %in% c("mu_b1", "mu_b3", "mu_b3_upr", "b4", "b5", "mu_d1") ~ "slope",
                            variable %in% c("mu_b2") ~ "slope_ele",
                            variable %in% c("mu_mid0") ~ "mid0", 
                            grepl("sigma", variable) ~ "sigma"),
           variable2 = variable,
           variable2 = case_when(variable == "b4" ~ "b3", 
                                 variable == "b5" ~ "b4", 
                                 variable == "mu_b3" ~ "mu_b5",
                                 variable == "mu_b3_upr" ~ "mu_b6", 
                                 variable == "sigma_b3" ~ "sigma_b5",
                                 variable == "sigma_b3_upr" ~ "sigma_b6", 
                                 variable == "mu_mid0" ~ "mu_mu",
                                 variable == "sigma_mid0" ~ "sigma_mu",
                                 TRUE ~ as.character(variable)),
           variable2 = factor(variable2, levels=lvls))


priors2 <- bind_rows(posterior, posterior2) %>%
    select(variable2, type) %>%
    unique %>% left_join(., priors) %>%
    mutate(variable_label = gsub("(.*)\\_(.*)", "\\2", variable2)) 

# plot
ggplot(posterior, aes(value)) + 
    facet_wrap(~variable2, scales="free", ncol=4) +
    scale_x_continuous(breaks = seq(-10, 10, 2), expand=c(0,0)) +
    geom_ribbon(data=priors2, aes(x, ymax=y), ymin=-Inf, fill="grey90", col="grey70") + 
    geom_density() + 
    geom_density(data=posterior2,lty=2) +
    theme_bw() +
    theme(strip.text = element_text(hjust=0, face="bold"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(),
          aspect.ratio = .6, 
          axis.text = element_text(colour="black")) +
    labs(x="Value", y = "Density")

ggsave("figures/posterior_prior_comps.png", units="cm", width=14, height=18)
