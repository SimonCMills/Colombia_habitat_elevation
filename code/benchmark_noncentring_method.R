# speed test of alternative non-centred parameterisation
# note: in initial runs I checked they were returning the same answer each time
# by binding along chains ..which all looked fine

library(posterior); library(dplyr); library(ggplot2); library(cmdstanr)

## simulate some data with a grouping factor
# note: n_pt is the number of observations within a group
n_pt <- 10
n_grp <- 30

# hyperparameters
mu_b0 <- 0
mu_b1 <- .5
sigma_b0 <- .5
sigma_b1 <- .5
rho <- .8
Rho <- matrix(c(1, rho, rho, 1), nrow=2)
sigmas <- c(sigma_b0, sigma_b1)
cov_ab <- sigma_b0 * sigma_b1 * rho
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)

# 3 versions (old_alt) doesn't store the transformed parameters, instead 
# defining them in the model block
mod_old <- cmdstan_model("code/stan_files/simple_2par_bernoulli.stan")
mod_old_alt <- cmdstan_model("code/stan_files/simple_2par_alt_bernoulli.stan")
mod_new <- cmdstan_model("code/stan_files/simple_2par_bernoulli_alternative_noncentring.stan")
                             

sim_fit <- replicate(100, simplify = F, {
    mvdraws <- MASS::mvrnorm(n_grp, c(mu_b0, mu_b1), Sigma)
    b0 <- mvdraws[,1]; b1 <- mvdraws[,2]
    
    dat <- data.frame(x = rnorm(n_pt*n_grp), 
                      id_grp = rep(1:n_grp, each=n_pt))
    dat2 <- within(dat, {
        mu <- boot::inv.logit(b0[id_grp] + b1[id_grp]*x)
        y <- rbinom(n_pt*n_grp, 1, mu)
    })
    
    # stan_data
    stan_data <- as.list(dat2)
    stan_data$n_pt <- nrow(dat2)
    stan_data$n_grp <- n_grp
    seed_i <- sample(1:1e7, 1)
    samples_old <- mod_old$sample(stan_data, seed=seed_i, chains=3, parallel_chains=3)
    samples_old_alt <- mod_old_alt$sample(stan_data, seed=seed_i, chains=3, parallel_chains=3)
    samples_new <- mod_new$sample(stan_data, seed=seed_i, chains=3, parallel_chains=3)
    
    tibble(t_old = mean(samples_old$time()$chains$total), 
           t_old_alt = mean(samples_old_alt$time()$chains$total), 
           t_new = mean(samples_new$time()$chains$total), 
           rel_diff = t_new/t_old)
})

# Save from previous run
# df_simfit <- bind_rows(sim_fit) 
# saveRDS(df_simfit, "data/speed_benchmark_1.rds")

    
ggplot(df_simfit, aes(log(t_new) - log(t_old))) + 
    geom_histogram(colour="black", fill="grey70") + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept = with(df_simfit, mean(log(t_new) - log(t_old))), lty=2) +
    theme_bw() +
    theme(aspect.ratio = .7, 
          axis.text = element_text(colour="black"), 
          panel.grid = element_blank()) +
    labs(title="Speed benchmark for noncentring; 50 reps, 30grp, 10pt per group (n=300 tot)")
ggsave("figures/speed_benchmark_1.png")

df_simfit2 <- bind_rows(sim_fit) %>% filter(complete.cases(.))
saveRDS(df_simfit2, "data/speed_benchmark_2.rds")
ggplot(df_simfit2, aes(log(t_new) - log(t_old))) + 
    geom_histogram(colour="black", fill="grey70") + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept = with(df_simfit, mean(log(t_new) - log(t_old))), lty=2) +
    theme_bw() +
    theme(aspect.ratio = .7, 
          axis.text = element_text(colour="black"), 
          panel.grid = element_blank()) +
    labs(title="Speed benchmark for noncentring; 100 reps, 30grp, 10pt per group (n=300 tot)", 
         caption = "Note: for some reason, final 3 reps failed failed by all three mods and were removed")

# ggsave("figures/speed_benchmark_1.png")

p1 <- ggplot(df_simfit2, aes(log(t_old_alt) - log(t_old))) + 
    geom_histogram(colour="black", fill="grey70") + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept = with(df_simfit2, mean(log(t_old_alt) - log(t_old))), lty=2) +
    theme_bw() +
    theme(#aspect.ratio = .7, 
          axis.text = element_text(colour="black"), 
          panel.grid = element_blank()) +
    labs(title="(a) Speed benchmark for old version, no-transformed-pars vs. old version", 
         caption = "Specification: 100 reps, 30grp, 10pt per group (n=300 tot) \n Note: for some reason, final 3 reps failed failed by all three mods and were removed", 
         x = "log(t), old-ver, no-transform pars  -  log(t), old version")

p2 <- ggplot(df_simfit2, aes(log(t_old_alt) - log(t_new))) + 
    geom_histogram(colour="black", fill="grey70") + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept = with(df_simfit2, mean(log(t_old_alt) - log(t_new))), lty=2) +
    theme_bw() +
    theme(#aspect.ratio = .7, 
          axis.text = element_text(colour="black"), 
          panel.grid = element_blank()) +
    labs(title="(b) Speed benchmark for old version, no-transformed-pars vs. new version", 
         caption = "Specification: 100 reps, 30grp, 10pt per group (n=300 tot) \n Note: for some reason, final 3 reps failed failed by all three mods and were removed", 
         x = "log(t), old-ver, no-transform pars  -  log(t), new version")

p3 <- ggplot(df_simfit2, aes(log(t_old) - log(t_new))) + 
    geom_histogram(colour="black", fill="grey70") + 
    geom_vline(xintercept=0) +
    geom_vline(xintercept = with(df_simfit2, mean(log(t_old) - log(t_new))), lty=2) +
    theme_bw() +
    theme(#aspect.ratio = .7, 
          axis.text = element_text(colour="black"), 
          panel.grid = element_blank()) +
    labs(title="(c) Speed benchmark for old version vs. new version", 
         caption = "Specification: 100 reps, 30grp, 10pt per group (n=300 tot) \n Note: for some reason, final 3 reps failed failed by all three mods and were removed", 
         x = "log(t), old-ver -  log(t), new version")

p_all <- egg::ggarrange(p1,p2,p3)
ggsave("figures/benchmarking_3_noncentring_methods.png", p_all, width=180, height=350, units="mm")



## confirm that these specifications are both doing the same thing
# library(posterior)
# draws_old_alt <- samples_old_alt$draws() %>%
#     subset_draws("raw|mu|sigma", regex=T) 
# 
# draws_old <- samples_old$draws() %>%
#     subset_draws("raw|mu|sigma", regex=T) 
# 
# bind_draws(draws_old_alt, draws_old, along="chain") %>% 
#     summarise_draws()
# 
# 
# draws_old <- samples_old$draws() %>%
#     subset_draws("b0\\[", regex=T)# %>% summarise_draws()
# 
# draws_new <- samples_new$draws() %>%
#     subset_draws("b0\\[", regex=T) # %>% summarise_draws()
# 
# bind_draws(draws_new, draws_old, along="chain") %>% 
#     summarise_draws()

