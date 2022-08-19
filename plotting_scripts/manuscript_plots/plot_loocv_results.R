library(rstan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)

TIME_CORRECTION_TYPE <<- 'beta'
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


plot_scatter = plot_observed_predicted_time(posterior_means, FALSE, TRUE, xlimits = c(0, 2.4), with_subject_legend = TRUE, with_fragment_legend = FALSE)
plot_scatter_loocv = plot_scatter + facet_grid(~fragment) + labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/loocv_scatter_by_frag.pdf'), plot = plot_scatter_loocv, width = 27, height = 9.5, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_loocv = plot_observed_predicted_time_histogram(posterior_means, FALSE, TRUE, xlimits = c(-2.1, 6))
ggsave(paste0('plots/manuscript_figs/loocv_hist_by_frag.pdf'), plot = plot_hist_loocv, width = 13, height = 11.25, units = 'in', dpi = 750, device = cairo_pdf)

