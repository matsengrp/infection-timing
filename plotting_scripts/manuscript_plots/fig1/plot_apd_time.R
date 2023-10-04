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

data[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
data[fragment != 1, fragment_long := paste0('Gene region ', fragment, ' (within pol)')]
data[infection_status == 'IN UTERO', infection_status_long := 'in-utero']
data[infection_status != 'IN UTERO', infection_status_long := 'after birth']

data[, count := .N, by = .(fragment, get(subject_var))]

data$infection_status_long = factor(data$infection_status_long, levels = c('in-utero', 'after birth'))

uteropalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(uteropalette) = unique(data[infection_status_long == 'in-utero'][[subject_var]])
afterpalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(afterpalette) = unique(data[infection_status_long == 'after birth'][[subject_var]])

palette = c(uteropalette, afterpalette)
palette_df = data.frame(palette)
palette_df[[subject_var]] = as.numeric(rownames(palette_df))

data = merge(data, palette_df, by = subject_var)

# Make lookup tables
subj_to_status = setNames(data$infection_status_long, as.factor(data[[subject_var]]))
status2shape = c("in-utero" = 16, "after birth" = 17)
status2line = c("in-utero" = 'solid', "after birth" = '41')
subj_to_palette = setNames(data$palette, as.factor(data[[subject_var]])) 

# # Make subj_to_status unique without dropping names
it_cols = c(subject_var, 'infection_status_long')
p_cols = c(subject_var, 'palette')
subj_to_status = subj_to_status[!duplicated(data[, ..it_cols])]
subj_to_status = subj_to_status[order(subj_to_status)]
subj_to_palette = subj_to_palette[!duplicated(data[, ..p_cols])]
subj_to_palette = subj_to_palette[order(subj_to_palette)]

valid = data[order(observed_time_since_infection)]

plot2 = ggplot(valid) +
    geom_point(aes(x = observed_time_since_infection, y = apd, color = as.factor(get(subject_var)), shape = infection_status_long), size = 5, alpha = 0.8) +
    geom_line(aes(x = observed_time_since_infection, y = apd, linetype = infection_status_long, color = as.factor(get(subject_var))), linewidth = 1.5, alpha = 0.8, key_glyph = 'abline') +
    facet_grid2(cols = vars(fragment_long), rows = vars(infection_status_long), scales = 'free_y', independent = 'y') +
    scale_color_manual(breaks = names(subj_to_status), guide  = guide_legend(ncol = 10, override.aes = list(linetype = status2line[subj_to_status], shape = status2shape[subj_to_status])), values = subj_to_palette[names(subj_to_status)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    scale_linetype_manual(values = status2line, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center", legend.key.width = unit(3.5, 'line')) +
    background_grid(major = 'xy') +
    xlab('\nTrue time since infection (months)\n') +
    ylab('Average pairwise diversity\n')+
    panel_border(color = 'gray60', size = 2)

name2 = paste0('plots/manuscript_figs/apd_time.pdf')
ggsave(name2, plot = plot2, width = 22, height = 14, units = 'in', dpi = 750, device = cairo_pdf)
