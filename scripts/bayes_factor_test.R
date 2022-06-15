library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(bridgesampling)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
NCPU <<- 2
PREPROCESS_DATA <<- TRUE
BAYES_FACTOR_VARIATION <<- args[2]
stopifnot(BAYES_FACTOR_VARIATION %in% c('no_fragment', 'no_subject', 'with_infection_time'))

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/bayes_factor_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)

# fit full model
model = fit_model(infant_data, total_iterations = 40000)

# fit model without indicated variation
NO_VAR_MODEL_FILE <<- str_replace(MODEL_FILE, '[.]stan', paste0('_', BAYES_FACTOR_VARIATION, '_variation.stan'))
no_var_model = fit_model(infant_data, model_file = NO_VAR_MODEL_FILE, total_iterations = 40000)

bayes_factor = compile_results(model, no_var_model)

