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
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))

infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
adult_estimate = fread('_ignore/neher_webtool_results.tsv')
subset = adult_estimate[, c('ptnum', 'fragment_int', 'year_visit', 'adult_model_time_since_infection_estimates')]
colnames(subset) = c('subject_id', 'fragment', 'observed_time_since_infection', 'adult_model_time_since_infection_estimates')

infant_data[, observed_time_since_infection := round(observed_time_since_infection, 6)]
subset[, observed_time_since_infection := round(observed_time_since_infection, 6)]
together = merge(infant_data, subset)

setnames(together, 'adult_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
plot_observed_predicted_time(together, FALSE, TRUE, FALSE, adult_data = TRUE)
plot_observed_predicted_time_histogram(together, FALSE, TRUE, adult_data = TRUE)
