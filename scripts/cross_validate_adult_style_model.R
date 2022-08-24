library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(L1pack)

TIME_CORRECTION_TYPE <<- 'none' 
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)

results = run_adult_style_model_loocv(infant_data)

file_name = get_adult_style_model_loocv_name()

fwrite(results, file_name, sep = '\t')
