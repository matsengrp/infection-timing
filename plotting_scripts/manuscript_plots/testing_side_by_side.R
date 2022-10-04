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
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

# get adult model results
name = str_split(TESTING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TESTING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'adult_model_estimates', paste0('adult_model_estimates_', name, '.csv'))

adult_results = fread(file_name)

true_results = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
true_results[, inftimeyears := inftimemonths/12]
true_results[, year_visit := month_visit/12]
true_results[, observed_time_since_infection := year_visit - inftimeyears]

adult_together = merge(adult_results, true_results, by= 'ptnum')
setnames(adult_together, 'adult_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
setnames(adult_together, 'fragment_int', 'fragment')
setnames(adult_together, 'ptnum', 'subject_id')
adult_together$model = 'adult model'
adult_mae = calculate_mae_dt(adult_together)
adult_mae$model = 'adult model'

# get adult-style model results
output_name = get_adult_style_model_predictions_name(TESTING_INFANT_DATA_PATH)
adult_style_results = fread(output_name)

adult_style_together = merge(adult_style_results, true_results, by.y= 'ptnum', by.x = 'subject_id')
setnames(adult_style_together, 'adult_style_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
adult_style_together$model = 'adult-style model'
adult_style_mae = calculate_mae_dt(adult_style_together)
adult_style_mae$model = 'adult-style model'

# get infant model results
file = get_model_output_filename(TESTING_INFANT_DATA_PATH)
posteriors = fread(file)
if ('mean_predicted_time_since_infection' %in% colnames(posteriors)){
    setnames(posteriors, 'mean_predicted_time_since_infection', 'mean_time_since_infection_estimate')
}

infant_together = merge(posteriors, true_results, by.y = 'ptnum', by.x = 'subject_id')
infant_together$model = 'infant hierarchical model'
infant_mae = calculate_mae_dt(infant_together)
infant_mae$model = 'infant hierarchical model'

together = rbind(adult_together, adult_style_together, infant_together, fill = TRUE)
mae = rbind(adult_mae, adult_style_mae, infant_mae) 

plot_scatter = plot_observed_predicted_time(together, FALSE, FALSE, write_plot = FALSE, with_subject_legend = TRUE, with_fragment_legend = FALSE)

plot_scatter_val = plot_scatter + facet_grid(cols = vars(fragment), rows = vars(model))+ labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/side_by_side_testing_scatter_by_frag.pdf'), plot = plot_scatter_val, width = 27, height = 27, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_val = plot_observed_predicted_time_histogram(together, FALSE, FALSE, write_plot = FALSE)+ 
    facet_grid(cols = vars(fragment), rows = vars(model))+
    geom_text(data = mae, x = 4, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 8)
ggsave(paste0('plots/manuscript_figs/side_by_side_testing_hist_by_frag.pdf'), plot = plot_hist_val, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

