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

output_name = get_adult_style_model_predictions_name(TESTING_INFANT_DATA_PATH)
results = fread(output_name)

true_results = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
true_results[, inftimeyears := inftimemonths/12]
true_results[, year_visit := month_visit/12]
true_results[, observed_time_since_infection := year_visit - inftimeyears]

together = merge(results, true_results, by.y= 'ptnum', by.x = 'subject_id')
setnames(together, 'adult_style_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')

mae = calculate_mae_dt(together)

plot_scatter = plot_observed_predicted_time(together, FALSE, FALSE, write_plot = FALSE, with_subject_legend = TRUE, with_fragment_legend = FALSE)

plot_scatter_val = plot_scatter + facet_grid(~fragment)+ labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/adult_style_model_testing_scatter_by_frag.pdf'), plot = plot_scatter_val, width = 27, height = 9.5, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_val = plot_observed_predicted_time_histogram(together, FALSE, FALSE, write_plot = FALSE, xlimits = c(-2, 2.25))+
    geom_text(data = mae, x = 2, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 8)
ggsave(paste0('plots/manuscript_figs/adult_style_model_testing_hist_by_frag.pdf'), plot = plot_hist_val, width = 13, height = 11.25, units = 'in', dpi = 750, device = cairo_pdf)

