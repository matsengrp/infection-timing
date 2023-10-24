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
results[, observed_time_since_infection := month_visit - inftimemonths]
setnames(results, 'adult_model_time_since_infection_estimates', 'mean')
setnames(results, 'fragment_int', 'fragment')
setnames(results, 'ptnum', 'subject_id')
results$model = 'adult-trained\nlinear models'
results$mean = results$mean *12
adult_mae = calculate_mae_dt(results, by_fragment = TRUE, var = 'mean')
adult_mae$model = 'adult-trained\nlinear models'

# get adult-style model results
file_name = get_adult_style_model_loocv_name()
adult_style_results = fread(file_name)
setnames(adult_style_results, 'adult_style_model_time_since_infection_estimates', 'mean')
adult_style_results$model = 'infant-trained\nlinear models'
adult_style_results$mean = adult_style_results$mean *12
adult_style_results$observed_time_since_infection = adult_style_results$observed_time_since_infection *12
adult_style_mae = calculate_mae_dt(adult_style_results, by_fragment = TRUE, var = 'mean')
adult_style_mae$model = 'infant-trained\nlinear models'

# get infant model results
name = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
file_name = file.path(OUTPUT_PATH, 'model_loocv', paste0('loocv_posteriors_', name, '_', TIME_CORRECTION_TYPE, '.csv'))
loocv = fread(file_name)
infant_data = as.data.table(configure_data(TRAINING_INFANT_DATA_PATH))
posterior_means = get_posterior_means(loocv, infant_data)
posterior_means$model = 'infant-trained\nhierarchical\nmodel'
setnames(posterior_means, 'mean_time_since_infection_estimate', 'mean')
posterior_means$mean = posterior_means$mean *12
posterior_means$observed_time_since_infection = posterior_means$observed_time_since_infection *12

infant_mae = calculate_mae_dt(posterior_means, by_fragment = TRUE, var = 'mean')
infant_mae$model = 'infant-trained\nhierarchical\nmodel'

together = rbind(results, adult_style_results, posterior_means, fill = TRUE)
mae = rbind(adult_mae, adult_style_mae, infant_mae)

together$difference = together$mean - together$observed_time_since_infection
together[fragment == 1, fragment_long := 'gene region 1 (within gag)']
together[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
mae[fragment == 1, fragment_long := 'gene region 1 (within gag)']
mae[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
together = merge(together, mae, by = c('fragment', 'fragment_long', 'model'))
plot_hist_val = ggplot(together) +
    geom_histogram(aes(x = difference), alpha = 0.7) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 2)+
    facet_grid(cols = vars(fragment_long), rows = vars(model))+
    geom_text(data = mae, x = 50, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12) +
    ylab('Observation count\n')+
    xlab('\nModel-derived time since infection - true time since infection (months)')

name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together_panelA.pdf')
ggsave(name, plot = plot_hist_val, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

cols = c('model', 'fragment_long', 'difference', 'mae')
plot_data = together[, ..cols]
colnames(plot_data) = c('model', 'gene_region', 'difference', 'mean_absolute_error')
name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together_panelA.csv')
fwrite(plot_data, name, sep = ',')

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)

intervals = get_posterior_prediction_coverage(loocv, actual_times = data)$data
intervals$mean = intervals[['q0.5']]
intervals$model = 'infant-trained\nhierarchical\nmodel'
intervals$mean = intervals$mean *12
intervals$observed_time_since_infection = intervals$observed_time_since_infection *12

adult_style_results$model = 'infant-trained\nlinear models'
tog = rbind(intervals, adult_style_results, fill = TRUE)

relation = data.table()
for (frag in c(1, 2, 3)){
    heir = lm(intervals[fragment == frag]$mean ~ intervals[fragment == frag]$observed_time_since_infection)
    lin = lm(adult_style_results[fragment == frag]$mean ~ adult_style_results[fragment == frag]$observed_time_since_infection)
    region = ' (within pol)'
    if (frag == 1){
        region = ' (within gag)'
    }
    temp_relation = data.table(model = c(paste0('infant-trained\nhierarchical\nmodel'), paste0('infant-trained\nlinear models')), r2 = c(round(summary(heir)$r.squared, 3), round(summary(lin)$r.squared, 3)), slope = c(round(coef(heir)[['intervals[fragment == frag]$observed_time_since_infection']], 3), round(coef(lin)[['adult_style_results[fragment == frag]$observed_time_since_infection']], 3)), intercept = c(round(coef(heir)[['(Intercept)']], 3), round(coef(lin)[['(Intercept)']], 3)))
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
    geom_text(data = relation, x = 6, y = 85, aes(label = paste0('R^2 = ', r2, '\nslope = ', slope, '\nintercept = ', intercept)), vjust = 2, size = 12) +
    geom_smooth(aes(x = observed_time_since_infection, y = mean), method = 'lm', color = 'gray60', size = 3, se = FALSE) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (months)') +
    ylab('Model-derived\ntime since infection (months)\n')+
    ylim(-12, 64)+
    panel_border(color = 'gray60', size = 2) 

name2 = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together_panelB.pdf')
ggsave(name2, plot = plot2, width = 20, height = 12, units = 'in', dpi = 750, device = cairo_pdf)
# saveRDS(plot2, paste0('plots/manuscript_figs/error_over_time.rds'))

cols = c('model', 'fragment_long', 'observed_time_since_infection', 'mean', 'q0.055', 'q0.945')
plot_data = tog[, ..cols]
colnames(plot_data) = c('model', 'gene_region', 'observed_time_since_infection', 'model_derived_time_since_infection_mean', 'model_derived_time_since_infection_0.055_quantile', 'model_derived_time_since_infection_0.945_quantile')
name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together_panelB.csv')
fwrite(plot_data, name, sep = ',')

tog[, resid := mean - observed_time_since_infection]
tog$model = str_replace_all(tog$model, '\n', ' ')
tog2 = tog[, mean(resid), by = .(model, observed_time_since_infection)]
heir2 = lm(tog2[model %like% 'hier']$V1 ~ tog2[model %like% 'hier']$observed_time_since_infection)
lin2 = lm(tog2[model %like% 'linear']$V1 ~ tog2[model %like% 'linear']$observed_time_since_infection)

plot3 = ggplot(tog2[order(observed_time_since_infection)], aes(x = observed_time_since_infection, y = V1, color = model, fill = model, group = model)) +
    geom_point(size = 7) +
    geom_smooth(method = 'lm', size = 3, se = FALSE) +
    # geom_line(aes(x = observed_time_since_infection, y = V1, color = model), size = 3) +
    geom_hline(yintercept = 0, size = 5, linetype = 'dashed', color = 'gray60') +
    theme_cowplot(font_family = 'Arial') +
    theme(legend.justification = 'center', legend.position = 'bottom', legend.direction = 'horizontal', axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 35), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (months)\n') +
    ylab('Mean residual\n(months)\n')+
    panel_border(color = 'gray60', size = 2) +
    scale_color_brewer(palette = 'Dark2') +
    scale_fill_brewer(palette = 'Dark2')

cols = c('model', 'observed_time_since_infection', 'V1')
plot_data = tog2[, ..cols]
colnames(plot_data) = c('model', 'observed_time_since_infection', 'model_residual')
name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together_panelC.csv')
fwrite(plot_data, name, sep = ',')


all = align_plots(plot2,plot_hist_val,plot3,  align = 'vh', axis = 'lr')

grid = plot_grid(all[[2]], NULL, all[[1]], NULL, all[[3]], nrow = 5, rel_heights = c(1.15, 0.05, 1.5, 0.05, 0.75), labels = c('A', '', 'B', '', 'C'), label_size = 40) 

name3 = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/fig5/training_error_together.pdf')
ggsave(name3, grid, width = 30, height = 38, units = 'in', dpi = 750, device = cairo_pdf)
# saveRDS(plot3, paste0('plots/manuscript_figs/residual.rds'))

print(summary(lm(intervals$mean ~ intervals$observed_time_since_infection)))
print(summary(lm(adult_style_results$mean ~ adult_style_results$observed_time_since_infection)))


