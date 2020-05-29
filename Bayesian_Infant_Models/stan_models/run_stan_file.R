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


fit1 <- stan(
  file = "multilevel_infant_model.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99), 
)

fit2 <- stan(
  file = "multilevel_infant_model_time_error.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99, max_treedepth = 20), 
)

save_models(fit2)

