library(rstan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)

args = commandArgs(trailingOnly = TRUE)

TIME_CORRECTION_TYPE <<- args[1]
stopifnot(TIME_CORRECTION_TYPE %in% c('uniform', 'beta', 'none'))
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))

infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
model = load_model_fit()

test_set_posteriors = predict_posterior(infant_data, model)
test_set_posterior_means = predict(test_set_posteriors) 

sim_data = simulate_apd_time_stan(model)

# plot just simulated data predicted time vs. apd
plot_apd_time(sim_data)

# plot simulated data and real data (observed time) with predicted time vs. apd
plot_apd_time(sim_data, infant_data)
