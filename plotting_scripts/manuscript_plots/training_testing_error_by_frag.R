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

###GET TESTING RESULTS
# get adult model results
name = str_split(TESTING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TESTING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'adult_model_estimates', paste0('adult_model_estimates_', name, '.csv'))

adult_results = fread(file_name)

true_results = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
true_results[, inftimeyears := inftimemonths/12]
true_results[, year_visit := month_visit/12]
true_results[, observed_time_since_infection := year_visit - inftimeyears]

test_adult_together = merge(adult_results, true_results, by= 'ptnum')
setnames(test_adult_together, 'adult_model_time_since_infection_estimates', 'mean')
setnames(test_adult_together, 'fragment_int', 'fragment')
setnames(test_adult_together, 'ptnum', 'subject_id')
test_adult_together$model = 'adult-trained\nlinear model'
test_adult_mae = calculate_mae_dt(test_adult_together, var = 'mean')
test_adult_mae$model = 'adult-trained\nlinear model'

# get adult-style model results
output_name = get_adult_style_model_predictions_name(TESTING_INFANT_DATA_PATH)
test_adult_style_results = fread(output_name)

test_adult_style_together = merge(test_adult_style_results, true_results, by.y= 'ptnum', by.x = 'subject_id')
setnames(test_adult_style_together, 'adult_style_model_time_since_infection_estimates', 'mean')
test_adult_style_together$model = 'infant-trained\nlinear model'
test_adult_style_mae = calculate_mae_dt(test_adult_style_together, var = 'mean')
test_adult_style_mae$model = 'infant-trained\nlinear model'

# get infant model results
file = get_model_output_filename(TESTING_INFANT_DATA_PATH)
posteriors = fread(file)
if ('mean_predicted_time_since_infection' %in% colnames(posteriors)){
    setnames(posteriors, 'mean_predicted_time_since_infection', 'mean')
}

test_infant_together = merge(posteriors, true_results, by.y = 'ptnum', by.x = 'subject_id')
test_infant_together$model = 'infant-trained\nhierarchical model'
setnames(test_infant_together, 'mean_time_since_infection_estimate', 'mean')

test_infant_mae = calculate_mae_dt(test_infant_together, var = 'mean')
test_infant_mae$model = 'infant-trained\nhierarchical model'

testing_together = rbind(test_adult_together, test_adult_style_together, test_infant_together, fill = TRUE)
testing_mae = rbind(test_adult_mae, test_adult_style_mae, test_infant_mae)
testing_together$type = 'testing data'
testing_mae$type = 'testing data'

#### get training results
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
setnames(results, 'adult_model_time_since_infection_estimates', 'mean')
setnames(results, 'fragment_int', 'fragment')
setnames(results, 'ptnum', 'subject_id')
results$model = 'adult-trained\nlinear model'
adult_mae = calculate_mae_dt(results, var = 'mean')
adult_mae$model = 'adult-trained\nlinear model'

# get adult-style model results
file_name = get_adult_style_model_loocv_name()
adult_style_results = fread(file_name)
setnames(adult_style_results, 'adult_style_model_time_since_infection_estimates', 'mean')
adult_style_results$model = 'infant-trained\nlinear model'
adult_style_mae = calculate_mae_dt(adult_style_results, var = 'mean')
adult_style_mae$model = 'infant-trained\nlinear model'

# get infant model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
loocv = fread(file_name)
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
posterior_means = get_posterior_means(loocv, infant_data)
posterior_means$model = 'infant-trained\nhierarchical model'
setnames(posterior_means, 'mean_time_since_infection_estimate', 'mean')
infant_mae = calculate_mae_dt(posterior_means, var = 'mean')
infant_mae$model = 'infant-trained\nhierarchical model'

training_together = rbind(results, adult_style_results, posterior_means, fill = TRUE)
training_mae = rbind(adult_mae, adult_style_mae, infant_mae)
training_together$type = 'training data'
training_mae$type = 'training data'

together = rbind(training_together, testing_together, fill = TRUE)
mae = rbind(training_mae, testing_mae, fill = TRUE)
together$difference = together$mean - together$observed_time_since_infection 
together[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
together[fragment == 2, fragment_long := 'Gene region 2 (within pol)']
together[fragment == 3, fragment_long := 'Gene region 3 (within pol)']
mae[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
mae[fragment == 2, fragment_long := 'Gene region 2 (within pol)']
mae[fragment == 3, fragment_long := 'Gene region 3 (within pol)']

together$type = factor(together$type, levels = c('training data', 'testing data'))
plot = ggplot(together, aes(x = difference, fill = type, group = type))+
    facet_grid(cols = vars(fragment_long), rows = vars(model))+
    geom_histogram(data = together, position = 'identity', alpha = 0.9) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal') +
    background_grid(major = 'xy') +
    geom_text(data = mae[type == 'testing data'], x = 18, y = 30, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12, color = '#d95f02') +
    geom_text(data = mae[type == 'training data'], x = 18, y = 35, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12, color = '#1b9e77') +
    ylab('Observation count\n')+
    xlab('ETI - TI (years)')+
    panel_border(color = 'gray60', size = 2)+
    scale_fill_brewer(palette = 'Dark2')+
    scale_color_brewer(palette = 'Dark2') +
    # guides(color = 'none') +
    labs(fill = 'Dataset')

ggsave(paste0('plots/manuscript_figs/side_by_side_hist_by_frag.pdf'), plot = plot, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)


