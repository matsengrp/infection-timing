library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(L1pack)

args = commandArgs(trailingOnly = TRUE)

DATA_PATH <<- args[1] 
TIME_CORRECTION_TYPE <<- 'none'
NCPU <<- 2
PREPROCESS_DATA <<- FALSE 

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

new_infant_data = as.data.table(data.frame(configure_newdata(DATA_PATH)))

results = data.table()
for (frag in unique(new_infant_data$fragment)){
    file_name = get_adult_style_model_name(frag)
    model = load_model_fit(file_name)

    data_subset = new_infant_data[fragment == frag]
    data_subset$adult_style_model_time_since_infection_estimates = stats::predict(model, newdata = data_subset, na.action = na.pass)
    data_subset$model_type = paste0(frag, '_fragment_model')
    results = rbind(results, data_subset)
}

output_name = get_adult_style_model_predictions_name(DATA_PATH)
# save posterior means
save_prediction_results(results, DATA_PATH, file_name = output_name)
