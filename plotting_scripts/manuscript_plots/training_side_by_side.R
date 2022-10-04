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
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

# get adult model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'adult_model_estimates', paste0('adult_model_estimates_', name, '.csv'))

results = fread(file_name)
results[, inftimeyears := inftimemonths/12]
results[, year_visit := month_visit/12]
results[, observed_time_since_infection := year_visit - inftimeyears]
setnames(results, 'adult_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
setnames(results, 'fragment_int', 'fragment')
setnames(results, 'ptnum', 'subject_id')
results$model = 'adult model'
adult_mae = calculate_mae_dt(results)
adult_mae$model = 'adult model'

# get adult-style model results
file_name = get_adult_style_model_loocv_name()
adult_style_results = fread(file_name)
setnames(adult_style_results, 'adult_style_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
adult_style_results$model = 'adult-style model'
adult_style_mae = calculate_mae_dt(adult_style_results)
adult_style_mae$model = 'adult-style model'

# get infant model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
loocv = fread(file_name)
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
posterior_means = get_posterior_means(loocv, infant_data)
posterior_means$model = 'infant hierarchical model'
infant_mae = calculate_mae_dt(posterior_means)
infant_mae$model = 'infant hierarchical model'

together = rbind(results, adult_style_results, posterior_means, fill = TRUE)
mae = rbind(adult_mae, adult_style_mae, infant_mae) 

plot_scatter = plot_observed_predicted_time(together, FALSE, FALSE, write_plot = FALSE, with_subject_legend = TRUE, with_fragment_legend = FALSE)

plot_scatter_val = plot_scatter + facet_grid(cols = vars(fragment), rows = vars(model))+ labs(color = 'ptnum') 
ggsave(paste0('plots/manuscript_figs/side_by_side_training_scatter_by_frag.pdf'), plot = plot_scatter_val, width = 27, height = 27, units = 'in', dpi = 750, device = cairo_pdf)

plot_hist_val = plot_observed_predicted_time_histogram(together, FALSE, FALSE, write_plot = FALSE)+ 
    facet_grid(cols = vars(fragment), rows = vars(model))+
    geom_text(data = mae, x = 4, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 8)
ggsave(paste0('plots/manuscript_figs/side_by_side_training_hist_by_frag.pdf'), plot = plot_hist_val, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

