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
true_results[, observed_time_since_infection := month_visit - inftimemonths]

adult_together = merge(adult_results, true_results, by= 'ptnum')
setnames(adult_together, 'adult_model_time_since_infection_estimates', 'mean')
adult_together$mean = adult_together$mean*12
setnames(adult_together, 'fragment_int', 'fragment')
setnames(adult_together, 'ptnum', 'subject_id')
adult_together$model = 'adult-trained\nlinear model'
adult_mae = calculate_mae_dt(adult_together, by_fragment = TRUE, var = 'mean')
adult_mae$model = 'adult-trained\nlinear model'

# get adult-style model results
output_name = get_adult_style_model_predictions_name(TESTING_INFANT_DATA_PATH)
adult_style_results = fread(output_name)

adult_style_together = merge(adult_style_results, true_results, by.y= 'ptnum', by.x = 'subject_id')
setnames(adult_style_together, 'adult_style_model_time_since_infection_estimates', 'mean')
adult_style_together$mean = adult_style_together$mean*12
adult_style_together$model = 'infant-trained\nlinear model'
adult_style_mae = calculate_mae_dt(adult_style_together, by_fragment = TRUE, var = 'mean')
adult_style_mae$model = 'infant-trained\nlinear model'

# get infant model results
# get infant model results
file = get_model_output_filename(TESTING_INFANT_DATA_PATH)
posteriors = fread(file)
if ('mean_predicted_time_since_infection' %in% colnames(posteriors)){
    setnames(posteriors, 'mean_predicted_time_since_infection', 'mean')
}

infant_together = merge(posteriors, true_results, by.y = 'ptnum', by.x = 'subject_id')
infant_together$model = 'infant-trained\nhierarchical model'
setnames(infant_together, 'mean_time_since_infection_estimate', 'mean')
infant_together$mean = infant_together$mean*12
infant_mae = calculate_mae_dt(infant_together, by_fragment = TRUE, var = 'mean')
infant_mae$model = 'infant-trained\nhierarchical model'

together = rbind(adult_together, adult_style_together, infant_together, fill = TRUE)
mae = rbind(adult_mae, adult_style_mae, infant_mae)

together$difference = together$mean - together$observed_time_since_infection 
together[fragment == 1, fragment_long := 'gene region 1 (within gag)']
together[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
mae[fragment == 1, fragment_long := 'gene region 1 (within gag)']
mae[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
together = merge(together, mae, by = c('fragment', 'fragment_long', 'model'))
plot = ggplot(together)+
    facet_grid(cols = vars(fragment_long), rows = vars(model))+
    geom_histogram(data = together, aes(x = difference), alpha = 0.7) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    geom_text(data = mae, x = 100, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12) +
    ylab('Observation count\n')+
    xlab('Model-derived time since infection - true time since infection (months)')+
    panel_border(color = 'gray60', size = 2)

name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig7/side_by_side_testing_hist_by_frag.pdf')
ggsave(name, plot = plot, width = 30, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

cols = c('model', 'fragment_long', 'difference', 'mae')
plot_data = together[, ..cols]
colnames(plot_data) = c('model', 'gene_region', 'difference', 'mean_absolute_error')
name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig7/side_by_side_testing_hist_by_frag.csv')
fwrite(plot_data, name, sep = ',')


