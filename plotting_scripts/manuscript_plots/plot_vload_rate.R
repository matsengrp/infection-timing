library(rstan)
library(rstanarm)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)

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

vload_data = configure_data(TRAINING_INFANT_DATA_PATH, include_vl = TRUE)
vload_data = as.data.table(vload_data)

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

tog = merge(vload_data, subset2)

tog$infection_status_long = factor(tog$infection_status_long, levels = c('in-utero', 'after birth'))

plot = ggplot(tog[order(observed_time_since_infection)])+
    geom_point(aes(y = vload, x = mean, color = as.factor(subject_id)), size = 5, alpha = 0.8) +
    facet_grid(rows = vars(fragment_long), cols = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center") +
    background_grid(major = 'xy') +
    ylab('log10(set-point viral load)\n') +
    xlab('\nMedian total APD slope\n')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=9))

name = paste0('plots/manuscript_figs/apd_slope_vload.pdf')
ggsave(name, plot = plot, width = 14, height = 20, units = 'in', dpi = 750, device = cairo_pdf)
