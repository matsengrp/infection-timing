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
setnames(adult_together, 'adult_model_time_since_infection_estimates', 'mean')
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
infant_mae = calculate_mae_dt(infant_together, by_fragment = TRUE, var = 'mean')
infant_mae$model = 'infant-trained\nhierarchical model'

together = rbind(adult_together, adult_style_together, infant_together, fill = TRUE)
mae = rbind(adult_mae, adult_style_mae, infant_mae)

together$difference = together$mean - together$observed_time_since_infection 
together[fragment == 1, fragment_long := 'gene region 1 (within gag)']
together[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]
mae[fragment == 1, fragment_long := 'gene region 1 (within gag)']
mae[fragment != 1, fragment_long := paste0('gene region ', fragment,' (within pol)')]

plot = ggplot()+
    facet_grid(cols = vars(fragment_long), rows = vars(model))+
    geom_histogram(data = together, aes(x = difference), alpha = 0.7) +
    geom_vline(xintercept = 0, color = 'black', size = 2) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    geom_text(data = mae, x = 17, y = Inf, aes(label = paste0('MAE = ', mae)), vjust = 2, size = 12) +
    ylab('Observation count\n')+
    xlab('Model-derived time since infection - true time since infection (years)')+
    panel_border(color = 'gray60', size = 2)

ggsave(paste0('plots/manuscript_figs/side_by_side_testing_hist_by_frag.pdf'), plot = plot, width = 30, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

# infant_data = configure_newdata(TESTING_INFANT_DATA_PATH)
# data = as.data.table(infant_data)
# actual = fread(TESTING_INFANT_DATA_TRUE_TIME_PATH)
# actual[, inftimeyears := inftimemonths/12]
# actual[, ptnum := as.character(ptnum)]

# model = load_model_fit()
# test_set_posterior_means = predict(infant_data, model) 
# subset2 = merge(test_set_posterior_means, actual, by.x = 'subject_id', by.y = 'ptnum')
# subset2[, observed_time_since_infection := ((month_visit)/12) - inftimeyears]
# setnames(subset2, '5.5%', 'cred_int_5')
# setnames(subset2, '94.5%', 'cred_int_95')
# subset2$model = 'infant-trained hierarchical model'

# adult_style_together$model = 'infant-trained linear model'
# adult_style_together[, subject_id := as.character(subject_id)]
# adult_style_together = merge(adult_style_together, actual, by.x = 'subject_id', by.y = 'ptnum')
# adult_style_together[, observed_time_since_infection := ((month_visit)/12) - inftimeyears]
# setnames(subset2, 'mean_predicted_time_since_infection', 'mean')
# tog = rbind(subset2, adult_style_together, fill = TRUE)
# tog[incat_hiv_bin == 'IN UTERO', infection_status_long := 'in-utero']
# tog[incat_hiv_bin != 'IN UTERO', infection_status_long := 'after birth']

# heir = lm(subset2$mean ~ subset2$observed_time_since_infection)
# lin = lm(adult_style_together$mean ~ adult_style_together$observed_time_since_infection)

# relation = data.table(model = c('infant-trained hierarchical model', 'infant-trained linear model'), r2 = c(round(summary(heir)$r.squared, 3), round(summary(lin)$r.squared, 3)), slope = c(round(coef(heir)[['subset2$observed_time_since_infection']], 3), round(coef(lin)[['adult_style_together$observed_time_since_infection']], 3)), intercept = c(round(coef(heir)[['(Intercept)']], 3), round(coef(lin)[['(Intercept)']], 3)))
# plot2 = ggplot(tog) +
#     facet_grid(cols = vars(model)) +
#     geom_abline(intercept = 0, size = 3, color = 'blue') +
#     geom_point(aes(x = observed_time_since_infection, y = mean), size = 7, alpha = 0.6) +
#     geom_segment(aes(x = observed_time_since_infection, y = cred_int_5, xend = observed_time_since_infection, yend = cred_int_95), size = 2, alpha = 0.5) +
#     geom_text(data = relation, x = 0.8, y = Inf, aes(label = paste0('R^2 = ', r2, '\nslope = ', slope, '\nintercept = ', intercept)), vjust = 2, size = 12) +
#     geom_smooth(aes(x = observed_time_since_infection, y = mean), method = 'lm', color = 'gray60', size = 3) +
#     theme_cowplot(font_family = 'Arial') +
#     theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
#     background_grid(major = 'xy') +
#     xlab('\nTI') +
#     ylab('ETI\n')+
#     panel_border(color = 'gray60', size = 2) 

# name2 = paste0('plots/manuscript_figs/testing_error_over_time.pdf')
# ggsave(name2, plot = plot2, width = 20, height = 12, units = 'in', dpi = 750, device = cairo_pdf)

# tog[, resid := mean - observed_time_since_infection]
# tog2 = tog[, mean(resid), by = .(model, observed_time_since_infection)]
# # heir2 = lm(tog2[model %like% 'hier']$V1 ~ tog2[model %like% 'hier']$observed_time_since_infection)
# # lin2 = lm(tog2[model %like% 'linear']$V1 ~ tog2[model %like% 'linear']$observed_time_since_infection)

# plot3 = ggplot(tog2[order(observed_time_since_infection)], aes(x = observed_time_since_infection, y = V1, color = model, fill = model, group = model)) +
#     geom_point(size = 7) +
#     geom_smooth(method = 'lm', size = 3) +
#     # geom_line(aes(x = observed_time_since_infection, y = V1, color = model), size = 3) +
#     geom_hline(yintercept = 0, size = 4, linetype = 'dashed', color = 'gray60') +
#     theme_cowplot(font_family = 'Arial') +
#     theme(legend.position = 'bottom', legend.direction = 'horizontal', axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 35), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
#     background_grid(major = 'xy') +
#     xlab('\nTI') +
#     ylab('Mean residual (ETI - TI)\n')+
#     panel_border(color = 'gray60', size = 2) +
#     scale_color_brewer(palette = 'Dark2') +
#     scale_fill_brewer(palette = 'Dark2')

# all = align_plots(plot2,plot,plot3,  align = 'vh', axis = 'lr')

# grid = plot_grid(all[[2]], NULL, all[[1]], NULL, all[[3]], nrow = 5, rel_heights = c(1.2, 0.1, 1, 0.1, 0.7), labels = c('A', '', 'B', '', 'C'), label_size = 35) 

# name3 = paste0('plots/manuscript_figs/testing_error_together.pdf')
# ggsave(name3, grid, width = 26, height = 36, units = 'in', dpi = 750, device = cairo_pdf)

# print(summary(lm(subset2$mean ~ subset2$observed_time_since_infection)))
# print(summary(lm(adult_style_together$mean ~ adult_style_together$observed_time_since_infection)))
