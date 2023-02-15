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
    ylab('APD slope (years/diversity)\n')+
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

data$infection_status_long = factor(data$infection_status_long, levels = c('in-utero', 'after birth'))

uteropalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(uteropalette) = unique(data[infection_status_long == 'in-utero']$subject_id)
afterpalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(afterpalette) = unique(data[infection_status_long == 'after birth']$subject_id)

palette = c(uteropalette, afterpalette)
palette_df = data.frame(palette)
palette_df$subject_id = as.numeric(rownames(palette_df))

data = merge(data, palette_df)

# Make lookup tables
subj_to_status = setNames(data$infection_status_long, as.factor(data$subject_id))
status2shape = c("in-utero" = 16, "after birth" = 17)
status2line = c("in-utero" = 'solid', "after birth" = '41')
subj_to_palette = setNames(data$palette, as.factor(data$subject_id)) 

# # Make subj_to_status unique without dropping names
subj_to_status = subj_to_status[!duplicated(data[, c("subject_id", "infection_status_long")])]
subj_to_status = subj_to_status[order(subj_to_status)]
subj_to_palette = subj_to_palette[!duplicated(data[, c("subject_id", "palette")])]
subj_to_palette = subj_to_palette[order(subj_to_palette)]

valid = data[order(observed_time_since_infection)]

plot2 = ggplot(valid) +
    geom_point(aes(x = observed_time_since_infection, y = apd, color = as.factor(subject_id), shape = infection_status_long), size = 5, alpha = 0.8) +
    geom_line(aes(x = observed_time_since_infection, y = apd, linetype = infection_status_long, color = as.factor(subject_id)), linewidth = 1.5, alpha = 0.8, key_glyph = 'abline') +
    facet_grid2(cols = vars(fragment_long), rows = vars(infection_status_long), scales = 'free_y', independent = 'y') +
    scale_color_manual(breaks = names(subj_to_status), guide  = guide_legend(ncol = 10, override.aes = list(linetype = status2line[subj_to_status], shape = status2shape[subj_to_status])), values = subj_to_palette[names(subj_to_status)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    scale_linetype_manual(values = status2line, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (years)\n') +
    ylab('Average pairwise diversity\n')+
    panel_border(color = 'gray60', size = 2)

name2 = paste0('plots/manuscript_figs/apd_time.pdf')
ggsave(name2, plot = plot2, width = 22, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

subj_to_status_temp = subj_to_status[names(subj_to_status) %in% as.character(unique(valid[fragment == 3]$subject_id))]
plot2_temp = ggplot(valid[fragment == 3]) +
    geom_point(aes(x = observed_time_since_infection, y = apd, color = as.factor(subject_id), shape = infection_status_long), size = 5, alpha = 0.8) +
    geom_line(aes(x = observed_time_since_infection, y = apd, linetype = infection_status_long, color = as.factor(subject_id)), linewidth = 1.5, alpha = 0.8, key_glyph = 'abline') +
    geom_abline(intercept = 0, slope = 0.008, color = 'black', linewidth = 3.5)+
    facet_grid2(cols = vars(infection_status_long)) +
    scale_color_manual(breaks = names(subj_to_status_temp), guide  = guide_legend(ncol = 6, override.aes = list(linetype = status2line[subj_to_status_temp], shape = status2shape[subj_to_status_temp])), values = subj_to_palette[names(subj_to_status_temp)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    scale_linetype_manual(values = status2line, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (years)\n') +
    ylab('Viral diversity\n')+
    panel_border(color = 'gray60', size = 2)

name2 = paste0('plots/manuscript_figs/apd_time_temp.pdf')
ggsave(name2, plot = plot2_temp, width = 15, height = 10.5, units = 'in', dpi = 750, device = cairo_pdf)

plot2_temp2 = ggplot(valid[fragment == 3]) +
    geom_point(aes(x = observed_time_since_infection, y = apd, color = as.factor(subject_id), shape = infection_status_long), size = 5, alpha = 0.8) +
    geom_line(aes(x = observed_time_since_infection, y = apd, linetype = infection_status_long, color = as.factor(subject_id)), linewidth = 1.5, alpha = 0.8, key_glyph = 'abline') +
    geom_abline(intercept = 0, slope = 0.008, color = 'black', linewidth = 3.5)+
    geom_abline(intercept = 0, slope = 0.002, color = 'red', linewidth = 3.5)+
    facet_grid2(cols = vars(infection_status_long)) +
    scale_color_manual(breaks = names(subj_to_status_temp), guide  = guide_legend(ncol = 6, override.aes = list(linetype = status2line[subj_to_status_temp], shape = status2shape[subj_to_status_temp])), values = subj_to_palette[names(subj_to_status_temp)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    scale_linetype_manual(values = status2line, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (years)\n') +
    ylab('Viral diversity\n')+
    panel_border(color = 'gray60', size = 2)

name2 = paste0('plots/manuscript_figs/apd_time_temp2.pdf')
ggsave(name2, plot = plot2_temp2, width = 15, height = 10.5, units = 'in', dpi = 750, device = cairo_pdf)



baseline = intervals[variable %like% 'baseline']

plot3 = ggplot()+
    geom_point(data = data[infection_status_long == 'in-utero'], aes(x = apd, y = observed_time_since_infection, color = as.factor(subject_id)), shape = 16, size = 5, alpha = 0.8) +
    geom_point(data = data[infection_status_long == 'after birth'], aes(x = apd, y = observed_time_since_infection, color = as.factor(subject_id)), shape = 17, size = 5, alpha = 0.8) +
    geom_abline(data = subset2, aes(intercept = 0, slope = mean, color = as.factor(subject_id), linetype = infection_status_long), linewidth = 1.5, alpha = 0.8)+
    facet_grid(cols = vars(fragment_long), rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    scale_color_manual(breaks = names(subj_to_status), guide  = guide_legend(ncol = 10, override.aes = list(linetype = status2line[subj_to_status], shape = status2shape[subj_to_status])), values = subj_to_palette[names(subj_to_status)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    scale_linetype_manual(values = status2line, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
    background_grid(major = 'xy') +
    ylab('True time since infection (years)\n') +
    xlab('\nAverage pairwise diversity\n')+
    panel_border(color = 'gray60', size = 2) 

name3 = paste0('plots/manuscript_figs/apd_time_model_slopes.pdf')
ggsave(name3, plot = plot3, width = 20, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

# baseline = intervals[variable %like% 'baseline']
# valid = valid[order(observed_time_since_infection)]

# cred_int = data.table(apd = seq(0, 0.0402, by = 0.0001))
# cred_int[, lower := baseline[['5.5%']] * apd]
# cred_int[, upper := baseline[['94.5%']] * apd]

# tog = data.table() 
# for (f in unique(valid$fragment_long)){
#     for (i in unique(valid$infection_status_long)){

#         temp = cred_int
#         temp$fragment_long = f

#         tog = rbind(tog, cred_int)

#     }
# }

# plot4 = ggplot(valid) +
#     geom_point(aes(y = observed_time_since_infection, x = apd, color = as.factor(subject_id), shape = infection_status_long), size = 5, alpha = 0.6) +
#     geom_path(aes(y = observed_time_since_infection, x = apd, linetype = infection_status_long, color = as.factor(subject_id)), linewidth = 1.5, alpha = 0.6, key_glyph = 'abline') +
#     geom_abline(data = baseline, aes(intercept = 0, slope = mean), color = 'black', linewidth = 2.5) +
#     geom_ribbon(data = cred_int, aes(x = apd, ymin = lower, ymax = upper), fill = 'gray', alpha = 0.3) +
#     facet_grid2(cols = vars(fragment_long), rows = vars(infection_status_long)) +
#     scale_color_manual(breaks = names(subj_to_status), guide  = guide_legend(ncol = 10, override.aes = list(linetype = status2line[subj_to_status], shape = status2shape[subj_to_status])), values = subj_to_palette[names(subj_to_status)], name = 'Individual') +
#     scale_shape_manual(values = status2shape, guide= "none")+
#     scale_linetype_manual(values = status2line, guide= "none")+
#     theme_cowplot(font_family = 'Arial') +
#     theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
#     background_grid(major = 'xy') +
#     ylab('\nTrue time since infection (years)\n') +
#     xlab('Average pairwise diversity\n')+
#     panel_border(color = 'gray60', size = 2)


# name4 = paste0('plots/manuscript_figs/apd_time_with_baseline.pdf')
# ggsave(name4, plot = plot4, width = 22, height = 14, units = 'in', dpi = 750, device = cairo_pdf)



