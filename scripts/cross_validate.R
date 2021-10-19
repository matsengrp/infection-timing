library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

options(mc.cores = NCPU)

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)

leave_one_out_posteriors = run_leave_one_out_validation(infant_data)

fwrite(leave_one_out_posteriors, paste0(OUTPUT_PATH, '/loocv_posteriors_', TIME_CORRECTION_TYPE, '.tsv'), sep = '\t')

# get posterior mean time since infection estimates
posterior_means = get_posterior_means(leave_one_out_posteriors, infant_data)

# get mean absolute difference
print(mean(posterior_means$time_abs_difference)
