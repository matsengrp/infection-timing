library(rstan)
library(rstanarm)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)
library(ggh4x)

TIME_CORRECTION_TYPE <<- 'beta'
NCPU <<- 2
PREPROCESS_DATA <<-  TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)

model = load_model_fit()
intervals = get_posterior_interval_data(model, prob = 0.89)

subset = intervals[variable %like% 'total_slope']
subset$number = sapply(subset$variable, function(x) str_split(x, '\\[')[[1]][2])
subset$number = as.numeric(str_remove(subset$number, '\\]'))
subset = subset[order(number)]

together = merge(subset, data, by = 'number')
cols = c('5.5%', '94.5%', 'mean', 'infection_status', 'fragment', 'subject_id')

subset2 = unique(together[, ..cols])
row_count = nrow(data[, .N, by = .(subject_id, fragment)])
stopifnot(nrow(subset2) == row_count)

subset2[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
subset2[fragment != 1, fragment_long := paste0('Gene region ', fragment, ' (within pol)')]
subset2[infection_status == 'IN UTERO', infection_status_long := 'in-utero']
subset2[infection_status != 'IN UTERO', infection_status_long := 'after birth']
setnames(subset2, '5.5%', 'cred_int_5')
setnames(subset2, '94.5%', 'cred_int_95')
subset2$converted_mean = 1/subset2$mean
subset2$converted_cred_int_5 = 1/subset2$cred_int_5
subset2$converted_cred_int_95 = 1/subset2$cred_int_95

# it = calculate_infection_time(TRAINING_INFANT_DATA_PATH, intervals)
# subject_order = it[order(it)]$subject_id
# subset2$subject_id = factor(subset2$subject_id, levels = subject_order)
subset2$subject_id = factor(subset2$subject_id)
subset2$infection_status_long = factor(subset2$infection_status_long, levels = c('in-utero', 'after birth'))
require(RColorBrewer)
it = unique(subset2$infection_status_long)
color_palette = c('in-utero' = '#1b9e77', 'after birth' = '#7570b3')

plot = ggplot(subset2) +
    # facet_grid2(rows = vars(fragment_long), cols = vars(infection_status_long), scales = "free_x", space = "free", independent = 'x')+
    facet_grid(rows = vars(fragment_long), cols = vars(infection_status_long), scales = 'free_x')+
    # facet_wrap(fragment_long ~ infection_status_long, scales="free_x") +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray60', size = 3) +
    geom_point(aes(x = subject_id, y = mean, color = infection_status_long), size = 7) +
    geom_segment(aes(x = subject_id, y = cred_int_5, xend = subject_id, yend = cred_int_95, color = infection_status_long), size = 2) +
    theme_cowplot(font_family = 'Arial') +
    theme(legend.position = 'none', axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5, hjust=1)) +
    background_grid(major = 'xy') +
    ylim(-10, 250) +
    xlab('\nIndividual') +
    ylab('Total APD slope (year/diversity)\n')+
    labs(color = 'Infection time') +
    panel_border(color = 'gray60', size = 2) +
    # scale_color_brewer(palette = 'Dark2') 
    scale_color_manual(values = color_palette)

name = paste0('plots/manuscript_figs/credible_intervals.pdf')
ggsave(name, plot = plot, width = 14, height = 20, units = 'in', dpi = 750, device = cairo_pdf)

data[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
data[fragment != 1, fragment_long := paste0('Gene region ', fragment, ' (within pol)')]
data[infection_status == 'IN UTERO', infection_status_long := 'in-utero']
data[infection_status != 'IN UTERO', infection_status_long := 'after birth']

data[, count := .N, by = .(fragment, subject_id)]

# subject_order = c(unique(data[infection_status_long == 'in-utero']$subject_id), unique(data[infection_status_long =='after birth']$subject_id))
# data$subject_id = factor(data$subject_id, levels = subject_order)
data$infection_status_long = factor(data$infection_status_long, levels = c('in-utero', 'after birth'))
plot2 = ggplot(data[count > 1][order(observed_time_since_infection)])+
    geom_point(aes(x = observed_time_since_infection, y = apd, color = as.factor(subject_id)), size = 5, alpha = 0.8) +
    # geom_smooth(aes(y = observed_time_since_infection, x = apd, color = as.factor(subject_id)), size = 1.5, alpha = 0.8, linetype = '21', method = 'lm', se = FALSE) +
    geom_line(aes(x = observed_time_since_infection, y = apd, color = as.factor(subject_id)), linewidth = 1.5, alpha = 0.8, linetype = '21') +
    facet_grid(cols = vars(fragment_long), rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3, 'line')) +
    background_grid(major = 'xy') +
    xlab('\nTime since infection (years)\n') +
    ylab('Average pairwise diversity\n')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=10))


name2 = paste0('plots/manuscript_figs/apd_time.pdf')
ggsave(name2, plot = plot2, width = 20, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

uteropalette = c(brewer.pal(n = 8, name = "Dark2"), "#1B9E77", "#D95F02", '#F781BF')
names(uteropalette) = unique(data[infection_status_long == 'in-utero']$subject_id)
afterpalette = c(brewer.pal(n = 8, name = "Dark2"), "#1B9E77", "#D95F02", '#F781BF')
names(afterpalette) = unique(data[infection_status_long == 'after birth']$subject_id)

palette = c(uteropalette, afterpalette)

plot3 = ggplot()+
    geom_point(data = data[count > 1][infection_status_long == 'in-utero'], aes(x = apd, y = observed_time_since_infection, color = as.factor(subject_id)), shape = 'circle', size = 5, alpha = 0.8) +
    geom_point(data = data[count > 1][infection_status_long == 'after birth'], aes(x = apd, y = observed_time_since_infection, color = as.factor(subject_id)), shape = 'triangle', size = 5, alpha = 0.8) +
    geom_abline(data = subset2[infection_status_long == 'in-utero'], aes(intercept = 0, slope = mean, color = as.factor(subject_id)), linewidth = 2, alpha = 0.8)+
    geom_abline(data = subset2[infection_status_long == 'after birth'], aes(intercept = 0, slope = mean, color = as.factor(subject_id)), linetype = '51', linewidth = 2, alpha = 0.8)+
    facet_grid(rows = vars(fragment_long), cols = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3, 'line')) +
    background_grid(major = 'xy') +
    ylab('Time since infection (years)\n') +
    xlab('\nAverage pairwise diversity\n')+
    scale_color_manual(values = palette) +
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=7))

name3 = paste0('plots/manuscript_figs/apd_time_model_slopes.pdf')
ggsave(name3, plot = plot3, width = 14, height = 20, units = 'in', dpi = 750, device = cairo_pdf)

# extra = intervals[variable %like% 'baseline' | variable %like% 'reparam']
# extra[, total_mean := mean]
# extra[variable %like% 'frag', total_mean := total_mean + extra[variable %like% 'baseline']$mean]
# plot4 = ggplot()+
#     geom_point(data = data[count > 1], aes(x = apd, y = observed_time_since_infection, color = as.factor(subject_id)), size = 5, alpha = 0.8) +
#     geom_abline(data = subset2, aes(intercept = 0, slope = mean, color = as.factor(subject_id)), linewidth = 2, alpha = 0.8)+
#     facet_grid(rows = vars(fragment_long), cols = vars(infection_status_long)) +
#     # facet_grid(cols = vars(fragment_long)) +
#     theme_cowplot(font_family = 'Arial') +
#     theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3, 'line')) +
#     background_grid(major = 'xy') +
#     ylab('Time since infection (years)\n') +
#     xlab('\nAverage pairwise diversity\n')+
#     labs(color = 'Individual') +
#     panel_border(color = 'gray60', size = 2) +
#     guides(color=guide_legend(ncol=7))



# print(intervals[variable %like% 'subject_slope_delta'])
# print(intervals[variable %like% 'fragment'])

# BAYES_FACTOR_VARIATION <<- 'with_infection_time'

# source('config/config.R')
# source(paste0(PROJECT_PATH, '/config/file_paths.R'))
# source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
# source(paste0(PROJECT_PATH, '/scripts/bayes_factor_functions.R'))

# # fit model without indicated variation
# NO_VAR_MODEL_FILE <<- str_replace(MODEL_FILE, '[.]stan', paste0('_', BAYES_FACTOR_VARIATION, '_variation.stan'))
# no_var_model = fit_model(infant_data, model_file = NO_VAR_MODEL_FILE, total_iterations = 40000)

# var_intervals = get_posterior_interval_data(no_var_model, prob = 0.89)
# print(var_intervals[variable %like% 'infection_time_slope_delta'])

all = align_plots(plot, plot2, align = 'vh', axis = 'tb')

grid = plot_grid(all[[2]], NULL, all[[1]], nrow = 3, rel_heights = c(1.2, 0.1, 1), labels = c('A', '', 'B'), label_size = 35) 

name3 = paste0('plots/manuscript_figs/apd_time_together.pdf')
ggsave(name3, grid, width = 20, height = 30, units = 'in', dpi = 750, device = cairo_pdf)

