library(rstan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)

TIME_CORRECTION_TYPE <<- 'beta'
NCPU <<- 2
PREPROCESS_DATA <<- FALSE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))

file = get_model_output_filename(TESTING_INFANT_DATA_PATH)
posteriors = fread(file)
if ('mean_predicted_time_since_infection' %in% colnames(posteriors)){
    setnames(posteriors, 'mean_predicted_time_since_infection', 'mean_time_since_infection_estimate')
}

true_results = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
true_results[, inftimeyears := inftimemonths/12]
true_results[, year_visit := month_visit/12]
true_results[, observed_time_since_infection := year_visit - inftimeyears]

together = merge(posteriors, true_results, by.y = 'ptnum', by.x = 'subject_id')
setnames(together, 'subject_id', 'subject')
plot = plot_apd_time_by_subject(together)
