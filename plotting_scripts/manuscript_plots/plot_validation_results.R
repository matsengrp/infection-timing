library(rstan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)

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
plot_scatter = plot_observed_predicted_time(together, FALSE, FALSE, write_plot = FALSE, xlimits = c(0, 2.4), ylimits = c(0, 7.4), with_subject_legend = TRUE, with_fragment_legend = FALSE)

plot_scatter_val = plot_scatter + facet_grid(~fragment)+ labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/testing_scatter_by_frag.pdf'), plot = plot_scatter_val, width = 27, height = 9.5, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_val = plot_observed_predicted_time_histogram(together, FALSE, FALSE, write_plot = FALSE, xlimits = c(-2.1, 6))
ggsave(paste0('plots/manuscript_figs/testing_hist_by_frag.pdf'), plot = plot_hist_val, width = 13, height = 11.25, units = 'in', dpi = 750, device = cairo_pdf)

