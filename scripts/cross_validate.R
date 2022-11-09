library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(doParallel)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none', 'beta_laplace'))
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

name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
dir.create(file.path(OUTPUT_PATH, 'model_loocv'))
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))


fwrite(leave_one_out_posteriors, file_name, sep = '\t')
