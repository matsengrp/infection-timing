library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(L1pack)

TIME_CORRECTION_TYPE <<- 'none' 
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/adult_style_model_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
infant_data_dt = as.data.table(data.frame(infant_data))

library(ggplot2)
library(cowplot)

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

new_infant_data = configure_newdata(TESTING_INFANT_DATA_PATH)
new_infant_data_dt = as.data.table(data.frame(new_infant_data))
true_results = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
true_results[, inftimeyears := inftimemonths/12]
true_results[, year_visit := month_visit/12]
true_results[, observed_time_since_infection := year_visit - inftimeyears]

testing = merge(new_infant_data_dt, true_results, by.x = 'subject_id', by.y = 'ptnum')
testing$type = 'testing data'
infant_data_dt$type = 'training data'

together = rbind(testing, infant_data_dt, fill = TRUE)
together[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
together[fragment == 2, fragment_long := 'Gene region 2 (within pol)']
together[fragment == 3, fragment_long := 'Gene region 3 (within pol)']

together$type = factor(together$type, levels = c('training data', 'testing data'))

plot_apd = ggplot(together, aes(x = apd, fill = type, group = type)) +
    facet_grid(rows = vars(fragment_long))+
    geom_histogram(data = together, position = 'identity', alpha = 0.9) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal') +
    background_grid(major = 'xy') +
    ylab('Observation count\n')+
    xlab('APD')+
    panel_border(color = 'gray60', size = 2)+
    scale_fill_brewer(palette = 'Dark2')+
    labs(fill = 'Dataset')

ggsave(paste0('plots/manuscript_figs/train_test_apd_hist.pdf'), plot = plot_apd, width = 18, height = 19, units = 'in', dpi = 750, device = cairo_pdf)

together_ti = unique(together[, c('subject_id', 'observed_time_since_infection', 'type')])
plot_ti = ggplot(together_ti, aes(x = observed_time_since_infection, fill = type, group = type)) +
    geom_histogram(data = together_ti, position = 'identity', alpha = 0.9) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position = 'bottom', legend.direction = 'horizontal') +
    background_grid(major = 'xy') +
    ylab('Observation count\n')+
    xlab('True time since infection (years)')+
    panel_border(color = 'gray60', size = 2)+
    scale_fill_brewer(palette = 'Dark2')+
    labs(fill = 'Dataset')

ggsave(paste0('plots/manuscript_figs/train_test_time_hist.pdf'), plot = plot_ti, width = 18, height = 6.5, units = 'in', dpi = 750, device = cairo_pdf)


