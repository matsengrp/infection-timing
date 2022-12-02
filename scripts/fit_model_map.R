library(data.table)
library(tidyverse)
library(rstan)
library(foreach)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
model = fit_model(infant_data, type = 'MAP')

save_model_fit(model, model_name = get_model_fit_name('MAP'))
