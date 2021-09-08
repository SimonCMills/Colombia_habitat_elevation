## run Stan model relating occupancy to elevation and habitat covariates. 
# read in relevant files, create zero-filled dataframe (i.e. add 0s where species
# are not observed on a given point-visit), run model & save.

# pass chain&thread info from bash
cpu_info <- as.numeric(commandArgs(T))
n_chains <- cpu_info[1]
n_threads <- cpu_info[2]

## Housekeeping ----
library("dplyr"); library("reshape2"); library("cmdstanr"); library("posterior");

# generate analysis dataset 
source("code/format_analysis_dataset.R")

# Run model ----
model_name <- "code/stan_files/model_elevation.stan"

# compile stan
mod <- cmdstan_model(model_name, 
                     cpp_options = list(stan_threads = T), 
                     force_recompile=TRUE)

# read previous fit for inverse metric and step size to speed adaptation
#prior_run <- read_cmdstan_csv("stan_temp/model_linear_elevation_effect_with_estimated_mid0_threads-202103091424-1-90e25d.csv")

# sample (writing to stan_temp)
samps <- mod$sample(data = stan_data, 
                    refresh = 50, 
                    chains = n_chains, 
                    parallel_chains = n_chains, 
                    #inv_metric = prior_run$inv_metric,
                    threads_per_chain = n_threads,
                    iter_warmup = 1000, 
                    iter_sampling = 1000,
                    output_dir = "stan_out", 
                    step_size = 0.02,
                    save_warmup = TRUE)

# copy files (renaming)
outputs <- samps$output_files()
model_name <- gsub("^(.*/.*)-2[0-9]+-.*", "\\1", outputs)
model_name_date <- paste0(model_name, "_", Sys.Date())
model_name_full <- paste0(model_name_date, "_chain", 1:length(model_name_date), ".csv")
file.copy(outputs, model_name_full)

# end
print("script end")
