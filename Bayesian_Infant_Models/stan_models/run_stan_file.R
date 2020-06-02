library(rstan)

setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/")
source("multilevel_infant_functions.R")
setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/stan_models")

data <- list(
    N = nrow(indexed_infant_data_cleaned_subset),
    N_subject_index = length(unique(indexed_infant_data_cleaned_subset$subject_index)),
    N_fragment_index = length(unique(indexed_infant_data_cleaned_subset$frag_index)),
    time = c(indexed_infant_data_cleaned_subset$time),
    apd = c(indexed_infant_data_cleaned_subset$apd), 
    subject_index = c(indexed_infant_data_cleaned_subset$subject_index), 
    fragment_index = c(indexed_infant_data_cleaned_subset$frag_index)
)


stan_varying_slopes_frag_subject <- stan(
  file = "multilevel_infant_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99), 
)

stan_time_correction_varying_slopes_frag_subject <- stan(
  file = "multilevel_infant_model_time_error.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99, max_treedepth = 20), 
)

stan_time_correction_varying_slopes_frag_subject_more_restricted <- stan(
  file = "multilevel_infant_model_time_error_more_restricted_cutoff.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99, max_treedepth = 20), 
)

stan_time_correction_non_linear <- stan(
  file = "multilevel_infant_model_time_error_non_linear.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99, max_treedepth = 20), 
)

setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/")
save_models("stan_time_correction_varying_slopes_frag_subject")

setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/")
save_models("stan_time_correction_varying_slopes_frag_subject_more_restricted")

