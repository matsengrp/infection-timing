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

cd4_data = configure_data(TRAINING_INFANT_DATA_PATH, include_meta_data = TRUE, meta_data_path = TRAINING_INFANT_METADATA_PATH)
cd4_data = as.data.table(cd4_data)

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

tog = merge(cd4_data, subset2)

tog$infection_status_long = factor(tog$infection_status_long, levels = c('in-utero', 'after birth'))

uteropalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(uteropalette) = unique(tog[infection_status_long == 'in-utero']$subject_id)
afterpalette = c(brewer.pal(n = 8, name = "Dark2"), "#E41A1C", "#377EB8", '#F781BF')
names(afterpalette) = unique(tog[infection_status_long == 'after birth']$subject_id)

palette = c(uteropalette, afterpalette)
palette = palette[!is.na(names(palette))]
palette_df = data.frame(palette)
palette_df$subject_id = as.numeric(rownames(palette_df))

tog = merge(tog, palette_df, by = 'subject_id')

# Make lookup tables
subj_to_status = setNames(tog$infection_status_long, as.factor(tog$subject_id))
status2shape = c("in-utero" = 16, "after birth" = 17)
subj_to_palette = setNames(tog$palette, as.factor(tog$subject_id)) 

# # Make subj_to_status unique without dropping names
subj_to_status = subj_to_status[!duplicated(tog[, c("subject_id", "infection_status_long")])]
subj_to_palette = subj_to_palette[!duplicated(tog[, c("subject_id", "palette")])]
subj_to_status = subj_to_status[order(subj_to_status)]
subj_to_palette = subj_to_palette[order(subj_to_palette)]

plot = ggplot(tog[order(observed_time_since_infection)])+
    geom_point(aes(x = avg_cd4, y = mean, shape = infection_status_long, color = as.factor(subject_id)), size = 5, alpha = 0.8) +
    facet_grid(cols = vars(fragment_long), rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    scale_color_manual(breaks = names(subj_to_status), guide  = guide_legend(ncol = 10, override.aes = list(shape = status2shape[subj_to_status])), values = subj_to_palette[names(subj_to_status)], name = 'Individual') +
    scale_shape_manual(values = status2shape, guide= "none")+
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center") +
    background_grid(major = 'xy') +
    xlab('\nRate of change in CD4+ T cell count\n') +
    ylab('Median APD slope (years/diversity)\n')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) 

name = paste0('plots/manuscript_figs/apd_slope_cd4percent.pdf')
ggsave(name, plot = plot, width = 20, height = 14, units = 'in', dpi = 750, device = cairo_pdf)
