library(data.table)
library(tidyverse)
library(rstan)
library(foreach)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
NCPU <<- 2
PREPROCESS_DATA <<-  FALSE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))

model = load_model_fit()

new_infant_data = configure_newdata(TESTING_INFANT_DATA_PATH)

test_set_posteriors = predict_posterior(new_infant_data, model)
test_set_posterior_means = predict(new_infant_data, model) 

# save posterior means
save_prediction_results(test_set_posterior_means, TESTING_INFANT_DATA_PATH)
