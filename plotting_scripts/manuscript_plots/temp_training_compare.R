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
results$model = 'adult-trained\nlinear model'
# get mean adult model results
adult_mean = results[, mean(mean_time_since_infection_estimate), by = .(year_visit, observed_time_since_infection, inftimeyears, incat_hiv, subject_id, model)]
setnames(adult_mean, 'V1', 'mean')
adult_mae = calculate_mae_dt(adult_mean, by_fragment = FALSE, var = 'mean')

# get adult-style model results
file_name = get_adult_style_model_loocv_name()
adult_style_results = fread(file_name)
setnames(adult_style_results, 'adult_style_model_time_since_infection_estimates', 'mean_time_since_infection_estimate')
adult_style_results$model = 'infant-trained\nlinear model'
# get mean adult model results
adult_style_mean = adult_style_results[, mean(mean_time_since_infection_estimate), by = .(observed_time_since_infection, subject_id, model)]
setnames(adult_style_mean, 'V1', 'mean')
adult_style_mae = calculate_mae_dt(adult_style_mean, by_fragment = FALSE, var = 'mean')

# get infant model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
loocv = fread(file_name)
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
posterior_means = get_posterior_means(loocv, infant_data)
posterior_means$model = 'infant-trained\nhierarchical model\n(MCMC)'
# get mean model results
infant_mean = posterior_means[, mean(mean_time_since_infection_estimate), by = .(observed_time_since_infection, subject_id, model)]
setnames(infant_mean, 'V1', 'mean')
infant_mae = calculate_mae_dt(infant_mean, by_fragment = FALSE, var = 'mean')

# get infant model results MAP
map_name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
map_name = str_replace(map_name, '.csv', '')
map_file_name = file.path(OUTPUT_PATH, 'model_loocv_map', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
map_loocv = fread(map_file_name)
map_posterior_means = get_posterior_means(map_loocv, infant_data)
map_posterior_means$model = 'infant-trained\nhierarchical model\n(MAP)'
# get mean model results
map_infant_mean = map_posterior_means[, mean(mean_time_since_infection_estimate), by = .(observed_time_since_infection, subject_id, model)]
setnames(map_infant_mean, 'V1', 'mean')
map_infant_mae = calculate_mae_dt(map_infant_mean, by_fragment = FALSE, var = 'mean')

together = rbind(adult_mean, adult_style_mean, infant_mean, map_infant_mean, fill = TRUE)
mae = data.table(mae = c(adult_mae, adult_style_mae, infant_mae, map_infant_mae), model = c('adult-trained\nlinear model', 'infant-trained\nlinear model', 'infant-trained\nhierarchical model\n(MCMC)', 'infant-trained\nhierarchical model\n(MAP)')) 

together$difference = together$mean - together$observed_time_since_infection

plot_hist_val = ggplot(together) +
    geom_histogram(aes(x = difference), alpha = 0.7) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 2)+
    facet_grid(rows = vars(model))+
    ylab('Observation count\n')+
    xlab('ETI - TI (years)')

plot_hist_val = plot_hist_val +
    geom_text(data = mae, x = 9, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12) 

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)

intervals = get_posterior_prediction_coverage(loocv, actual_times = data)$data
intervals$mean = intervals[['q0.5']]
intervals$model = 'infant-trained hierarchical model (MCMC)'
adult_style_mean$model = 'infant-trained linear model'
map_infant_mean$model = 'infant-trained hierarchical model (MAP)'
tog = rbind(map_infant_mean, intervals, adult_style_mean, fill = TRUE)

heir = lm(intervals$mean ~ intervals$observed_time_since_infection)
heir2 = lm(map_infant_mean$mean ~ map_infant_mean$observed_time_since_infection)
lin = lm(adult_style_mean$mean ~ adult_style_mean$observed_time_since_infection)

relation = data.table(model = c('infant-trained hierarchical model (MCMC)', 'infant-trained hierarchical model (MAP)', 'infant-trained linear model'), r2 = c(round(summary(heir)$r.squared, 3), round(summary(heir2)$r.squared, 3),round(summary(lin)$r.squared, 3)), slope = c(round(coef(heir)[['intervals$observed_time_since_infection']], 3), round(coef(heir2)[['map_infant_mean$observed_time_since_infection']], 3), round(coef(lin)[['adult_style_mean$observed_time_since_infection']], 3)), intercept = c(round(coef(heir)[['(Intercept)']], 3), round(coef(heir2)[['(Intercept)']], 3),round(coef(lin)[['(Intercept)']], 3)))
plot2 = ggplot(tog) +
    facet_grid(cols = vars(model)) +
    geom_abline(intercept = 0, size = 3, color = 'blue') +
    geom_point(aes(x = observed_time_since_infection, y = mean), size = 7, alpha = 0.6) +
    geom_segment(aes(x = observed_time_since_infection, y = q0.05, xend = observed_time_since_infection, yend = q0.95), size = 2, alpha = 0.5) +
    geom_text(data = relation, x = 0.4, y = Inf, aes(label = paste0('R^2 = ', r2, '\nslope = ', slope, '\nintercept = ', intercept)), vjust = 2, size = 12) +
    geom_smooth(aes(x = observed_time_since_infection, y = mean), method = 'lm', color = 'gray60', size = 3) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    xlab('\nTI') +
    ylab('ETI\n')+
    panel_border(color = 'gray60', size = 2) 

all = align_plots(plot2,plot_hist_val,align = 'vh', axis = 'lr')

grid = plot_grid(all[[2]], NULL, all[[1]], nrow = 3, rel_heights = c(1.2, 0.05, 0.8), labels = c('A', '', 'B'), label_size = 35) 

name3 = paste0('plots/manuscript_figs/training_error_together_with_map.pdf')
ggsave(name3, grid, width = 32, height = 25, units = 'in', dpi = 750, device = cairo_pdf)

print(summary(lm(intervals$mean ~ intervals$observed_time_since_infection)))
print(summary(lm(adult_style_results$mean ~ adult_style_results$observed_time_since_infection)))
