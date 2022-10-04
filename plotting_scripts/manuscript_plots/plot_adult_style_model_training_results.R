library(rstan)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)

TIME_CORRECTION_TYPE <<- 'none'
NCPU <<- 2
PREPROCESS_DATA <<- FALSE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

file_name = get_adult_style_model_loocv_name()
results = fread(file_name)

setnames(results, 'adult_style_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')

mae = calculate_mae_dt(results)

plot_scatter = plot_observed_predicted_time(results, FALSE, FALSE, write_plot = FALSE, with_subject_legend = TRUE, with_fragment_legend = FALSE)

plot_scatter_val = plot_scatter + facet_grid(~fragment)+ labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/adult_style_model_loocv_scatter_by_frag.pdf'), plot = plot_scatter_val, width = 27, height = 9.5, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_val = plot_observed_predicted_time_histogram(results, FALSE, FALSE, write_plot = FALSE, xlimits = c(-2, 2.25))+
    geom_text(data = mae, x = 2, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 8)
ggsave(paste0('plots/manuscript_figs/adult_style_model_loocv_hist_by_frag.pdf'), plot = plot_hist_val, width = 13, height = 11.25, units = 'in', dpi = 750, device = cairo_pdf)

