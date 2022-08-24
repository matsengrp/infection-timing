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
infant_data_dt = as.data.table(data.frame(infant_data))

for (frag in unique(infant_data$fragment)){
    data_subset = infant_data_dt[fragment == frag]
    subset_model = lad(observed_time_since_infection ~ apd, data = data_subset, na.action = na.pass)

    file_name = get_adult_style_model_name(frag)

    save_model_fit(subset_model, file_name)
} 
