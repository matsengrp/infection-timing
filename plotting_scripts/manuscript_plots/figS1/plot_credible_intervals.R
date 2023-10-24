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
MAP_SUBJECTS <<- TRUE

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
if (isTRUE(MAP_SUBJECTS)){
    map = fread('_ignore/subject_mapping.tsv')
    data = merge(data, map, by = 'subject_id')
    subject_var = 'subject_map'
} else {
    subject_var = 'subject_id'
}
data[, observed_time_since_infection := observed_time_since_infection * 12]

model = load_model_fit()
intervals = get_posterior_interval_data(model, prob = 0.89)

subset = intervals[variable %like% 'total_slope']
subset$number = sapply(subset$variable, function(x) str_split(x, '\\[')[[1]][2])
subset$number = as.numeric(str_remove(subset$number, '\\]'))
subset = subset[order(number)]

together = merge(subset, data, by = 'number')
cols = c('5.5%', '94.5%', 'mean', 'infection_status', 'fragment', subject_var)

subset2 = unique(together[, ..cols])
row_count = nrow(data[, .N, by = .(get(subject_var), fragment)])
stopifnot(nrow(subset2) == row_count)

subset2[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
subset2[fragment != 1, fragment_long := paste0('Gene region ', fragment, ' (within pol)')]
subset2[infection_status == 'IN UTERO', infection_status_long := 'in-utero']
subset2[infection_status != 'IN UTERO', infection_status_long := 'after birth']
setnames(subset2, '5.5%', 'cred_int_5')
setnames(subset2, '94.5%', 'cred_int_95')
subset2[, mean := mean*12]
subset2[, cred_int_5 := cred_int_5*12]
subset2[, cred_int_95 := cred_int_95*12]

subset2$converted_mean = 1/subset2$mean
subset2$converted_cred_int_5 = 1/subset2$cred_int_5
subset2$converted_cred_int_95 = 1/subset2$cred_int_95

subset2[[subject_var]] = factor(subset2[[subject_var]])
subset2$infection_status_long = factor(subset2$infection_status_long, levels = c('in-utero', 'after birth'))
require(RColorBrewer)
it = unique(subset2$infection_status_long)
color_palette = c('in-utero' = '#1b9e77', 'after birth' = '#7570b3')

plot = ggplot(subset2) +
    # facet_grid2(rows = vars(fragment_long), cols = vars(infection_status_long), scales = "free_x", space = "free", independent = 'x')+
    facet_grid(rows = vars(fragment_long), cols = vars(infection_status_long), scales = 'free_x')+
    # facet_wrap(fragment_long ~ infection_status_long, scales="free_x") +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray60', size = 3) +
    geom_point(aes(x = get(subject_var), y = mean, color = infection_status_long), size = 7) +
    geom_segment(aes(x = get(subject_var), y = cred_int_5, xend = get(subject_var), yend = cred_int_95, color = infection_status_long), size = 2) +
    theme_cowplot(font_family = 'Arial') +
    theme(legend.position = 'none', axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5, hjust=1)) +
    background_grid(major = 'xy') +
    ylim(-120, 3000) +
    xlab('\nIndividual') +
    ylab('APD slope (months/diversity)\n')+
    labs(color = 'Infection time') +
    panel_border(color = 'gray60', size = 2) +
    # scale_color_brewer(palette = 'Dark2') 
    scale_color_manual(values = color_palette)

name = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/figS1/credible_intervals.pdf')
ggsave(name, plot = plot, width = 14, height = 20, units = 'in', dpi = 750, device = cairo_pdf)

cols = c(subject_var, 'infection_status_long', 'fragment_long', 'mean', 'cred_int_5', 'cred_int_95')
plot_data = subset2[, ..cols]
colnames(plot_data) = c('individual', 'mode_of_infection', 'gene_region', 'APD_slope_mean', 'APD_slope_0.055_quantile', 'APD_slope_0.945_quantile')
name2 = paste0(PROJECT_PATH, '/plotting_scripts/manuscript_plots/figS1/credible_intervals.csv')

fwrite(plot_data, name2, sep = ',')

