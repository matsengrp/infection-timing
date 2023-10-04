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

# get infant model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
loocv = fread(file_name)
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
posterior_means = get_posterior_means(loocv, infant_data)
setnames(posterior_means, 'mean_time_since_infection_estimate', 'mean')
posterior_means$mean = posterior_means$mean*12
posterior_means$model = 'infant-trained\nhierarchical\nmodel (Normal)'
infant_mae = calculate_mae_dt(posterior_means, by_fragment = TRUE, var = 'mean')
infant_mae$model = 'infant-trained\nhierarchical\nmodel (Normal)'

# get infant model results 
Laplace_name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
Laplace_name = str_replace(Laplace_name, '.csv', '')
Laplace_file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '_laplace.csv'))
Laplace_loocv = fread(Laplace_file_name)
Laplace_posterior_means = get_posterior_means(Laplace_loocv, infant_data)
setnames(Laplace_posterior_means, 'mean_time_since_infection_estimate', 'mean')
Laplace_posterior_means$mean = Laplace_posterior_means$mean*12
Laplace_posterior_means$model = 'infant-trained\nhierarchical\nmodel (Laplace)'
Laplace_infant_mae = calculate_mae_dt(Laplace_posterior_means, by_fragment = TRUE, var = 'mean')
Laplace_infant_mae$model = 'infant-trained\nhierarchical\nmodel (Laplace)'

together = rbind(posterior_means, Laplace_posterior_means, fill = TRUE)
together$observed_time_since_infection = together$observed_time_since_infection*12
mae = rbind(infant_mae, Laplace_infant_mae)

together$difference = together$mean - together$observed_time_since_infection
together[fragment == 1, fragment_long := 'gene region 1 (within gag)']
together[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
mae[fragment == 1, fragment_long := 'gene region 1 (within gag)']
mae[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]


plot_hist_val = ggplot(together) +
    geom_histogram(aes(x = difference), alpha = 0.7) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 2)+
    facet_grid(cols = vars(fragment_long), rows = vars(model))+
    ylab('Observation count\n')+
    xlab('Model-derived time since infection - true time since infection (months)')

plot_hist_val = plot_hist_val +
    geom_text(data = mae, x = 12, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12) 

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)
data$observed_time_since_infection = data$observed_time_since_infection*12

intervals = get_posterior_prediction_coverage(loocv, actual_times = data)$data
intervals$mean = intervals[['q0.5']]
intervals$mean = intervals$mean *12

intervals$model = 'infant-trained\nhierarchical\nmodel (Normal)'
Laplace_intervals = get_posterior_prediction_coverage(Laplace_loocv, actual_times = data)$data
Laplace_intervals$mean = Laplace_intervals[['q0.5']]
Laplace_intervals$mean = Laplace_intervals$mean *12
Laplace_intervals$model = 'infant-trained\nhierarchical\nmodel (Laplace)'

tog = rbind(Laplace_intervals, intervals, fill = TRUE)

relation = data.table()
for (frag in c(1, 2, 3)){
    heir = lm(intervals[fragment == frag]$mean ~ intervals[fragment == frag]$observed_time_since_infection)
    heir2 = lm(Laplace_intervals[fragment == frag]$mean ~ Laplace_intervals[fragment == frag]$observed_time_since_infection)

    region = ' (within pol)'
    if (frag == 1){
        region = ' (within gag)'
    }
    temp_relation = data.table(model = c(paste0('infant-trained\nhierarchical\nmodel (Normal)'), paste0('infant-trained\nhierarchical\nmodel (Laplace)')), r2 = c(round(summary(heir)$r.squared, 3), round(summary(heir2)$r.squared, 3)), slope = c(round(coef(heir)[['intervals[fragment == frag]$observed_time_since_infection']], 3), round(coef(heir2)[['Laplace_intervals[fragment == frag]$observed_time_since_infection']], 3)), intercept = c(round(coef(heir)[['(Intercept)']], 3), round(coef(heir2)[['(Intercept)']], 3)))
    temp_relation[, fragment_long := paste0('gene region ', frag, region)]
    relation = rbind(temp_relation, relation)
    tog[fragment == frag, fragment_long := paste0('gene region ', frag, region)]
}

tog[q0.945 > 5.5, q0.945 := Inf]


plot2 = ggplot(tog) +
    facet_grid(rows = vars(model), cols = vars(fragment_long)) +
    geom_abline(intercept = 0, size = 3, color = 'blue') +
    geom_point(aes(x = observed_time_since_infection, y = mean), size = 7, alpha = 0.6) +
    geom_segment(aes(x = observed_time_since_infection, y = q0.055*12, xend = observed_time_since_infection, yend = q0.945*12), size = 1.5, alpha = 0.3) +
    geom_text(data = relation, x = 5, y = 52, aes(label = paste0('R^2 = ', r2, '\nslope = ', slope, '\nintercept = ', intercept)), size = 12) +
    geom_smooth(aes(x = observed_time_since_infection, y = mean), method = 'lm', color = 'gray60', size = 3) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (months)') +
    ylab('Model-derived\ntime since infection (months)\n')+
    panel_border(color = 'gray60', size = 2) 

all = align_plots(plot2,plot_hist_val,align = 'vh', axis = 'lr')

grid = plot_grid(all[[2]], NULL, all[[1]], nrow = 3, rel_heights = c(1, 0.05, 1.6), labels = c('A', '', 'B'), label_size = 35) 

name3 = paste0('plots/manuscript_figs/training_error_together_with_Laplace.pdf')
ggsave(name3, grid, width = 40, height = 26, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

print(summary(lm(intervals$mean ~ intervals$observed_time_since_infection)))
