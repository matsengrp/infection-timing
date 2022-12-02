library(data.table)
library(tidyverse)
library(rstan)
library(foreach)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
DATA_PATH <<- args[2] 
NCPU <<- 2
PREPROCESS_DATA <<- FALSE 

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))

model = load_model_fit(model_name = get_model_fit_name('MAP'))

new_infant_data = configure_newdata(DATA_PATH)

test_set_posterior_means = predict(new_infant_data, model, type = 'MAP') 

# save posterior means
save_prediction_results(test_set_posterior_means, DATA_PATH, type = 'MAP')
