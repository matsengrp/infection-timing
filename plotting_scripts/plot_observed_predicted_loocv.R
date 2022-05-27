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

loocv = fread(paste0(OUTPUT_PATH, '/loocv_posteriors_', TIME_CORRECTION_TYPE, '.tsv'))
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))

# get posterior mean time since infection estimates
posterior_means = get_posterior_means(loocv, infant_data)


plot_observed_predicted_time(posterior_means, FALSE, TRUE, FALSE)
plot_observed_predicted_time_histogram(posterior_means, FALSE, TRUE)
