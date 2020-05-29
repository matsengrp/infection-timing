library(rstan)

source("multilevel_infant_functions.R")
#source("time_error_infant_model_scratch.stan")

data <- list(
    N = nrow(indexed_infant_data_cleaned_subset),
    time = c(indexed_infant_data_cleaned_subset$time),
    apd = c(indexed_infant_data_cleaned_subset$apd)
)


fit1 <- stan(
  file = "time_error_infant_model_scratch.stan",  # Stan program
  data = data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 10000,          # number of warmup iterations per chain
  iter = 20000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  control = list(adapt_delta = 0.99), 
)