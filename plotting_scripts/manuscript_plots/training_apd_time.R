library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(cowplot)

TIME_CORRECTION_TYPE <<- 'beta'
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)

data[fragment == 1, fragment_long := 'Gene region 1 (within gag)']
data[fragment != 1, fragment_long := paste0('Gene region ', fragment, ' (within pol)')]
data[infection_status == 'IN UTERO', infection_status_long := 'in-utero']
data[infection_status != 'IN UTERO', infection_status_long := 'after birth']

subject_order = c(unique(data[infection_status_long == 'in-utero']$subject_id), unique(data[infection_status_long =='after birth']$subject_id))
data$subject_id = factor(data$subject_id, levels = subject_order)

plot = ggplot(data[order(observed_time_since_infection)])+
    geom_point(aes(x = observed_time_since_infection, y = apd, color = subject_id), size = 5, alpha = 0.8) +
    # geom_line(aes(x = observed_time_since_infection, y = apd, color = subject_id), size = 1.5, alpha = 0.8, linetype = 'dashed') +
    geom_line(aes(x = observed_time_since_infection, y = apd, color = subject_id), size = 1.5, alpha = 0.8, linetype = '21') +
    facet_grid(cols = vars(fragment_long), rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center") +
    background_grid(major = 'xy') +
    xlab('TI\n') +
    ylab('APD')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=11))

name = paste0('plots/manuscript_figs/apd_time_raw.pdf')
ggsave(name, plot = plot, width = 20, height = 18, units = 'in', dpi = 750, device = cairo_pdf)

# get VL data
all_data = fread(TRAINING_INFANT_DATA_PATH)
cols = c('ptnum', 'month_visit', 'inftimemonths', 'vload', 'incat_hiv')
all_data = unique(all_data[, ..cols])
all_data = all_data[vload != 'NA' & vload != -99]
all_data[, observed_time_since_infection := month_visit/12]
all_data[, inftimeyears := inftimemonths/12]
all_data[incat_hiv != 'IN UTERO', observed_time_since_infection := observed_time_since_infection - inftimeyears]

all_data$ptnum = as.factor(all_data$ptnum)
tog = merge(data, all_data, by.x = c('subject_id', 'observed_time_since_infection'), by.y = c('ptnum', 'observed_time_since_infection'), all.x = TRUE)
# get coefs of determination

plot2 = ggplot(tog[!is.na(vload) & !(vload == -99)][order(observed_time_since_infection)])+
    geom_point(aes(x = observed_time_since_infection, y = log10(vload), color = subject_id), size = 5, alpha = 0.8) +
    # geom_line(aes(x = observed_time_since_infection, y = vload, color = subject_id), size = 1.5, alpha = 0.8, linetype = 'dashed') +
    geom_line(aes(x = observed_time_since_infection, y = log10(vload), color = subject_id), size = 1.5, alpha = 0.8, linetype = '21') +
    facet_grid(rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center") +
    background_grid(major = 'xy') +
    xlab('TI\n') +
    ylab('log10(Viral Load)')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=11))

name2 = paste0('plots/manuscript_figs/vl_time_raw.pdf')
ggsave(name2, plot = plot2, width = 8, height = 18, units = 'in', dpi = 750, device = cairo_pdf)

plot3 = ggplot(tog[!is.na(vload) & !(vload == -99)][order(observed_time_since_infection)])+
    geom_point(aes(x = apd, y = log10(vload), color = subject_id), size = 5, alpha = 0.8) +
    # geom_line(aes(x = observed_time_since_infection, y = vload, color = subject_id), size = 1.5, alpha = 0.8, linetype = 'dashed') +
    geom_line(aes(x = apd, y = log10(vload), color = subject_id), size = 1.5, alpha = 0.8, linetype = '21') +
    facet_grid(cols = vars(fragment_long), rows = vars(infection_status_long)) +
    # facet_grid(cols = vars(fragment_long)) +
    theme_cowplot(font_family = 'Arial') +
    theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.position="bottom", legend.direction="horizontal", legend.justification="center") +
    background_grid(major = 'xy') +
    xlab('APD\n') +
    ylab('log10(Viral Load)')+
    labs(color = 'Individual') +
    panel_border(color = 'gray60', size = 2) +
    guides(color=guide_legend(ncol=11))

name3 = paste0('plots/manuscript_figs/vl_apd_raw.pdf')
ggsave(name3, plot = plot3, width = 20, height = 18, units = 'in', dpi = 750, device = cairo_pdf)


r = data.table()
data[fragment == 1, fragment_long := 'Gene region 1\n(within gag)']
data[fragment != 1, fragment_long := paste0('Gene region ', fragment, '\n(within pol)')]

for (sub in unique(data$subject_id)){
    for (frag in unique(data$fragment_long)){
        temp = data[subject_id == sub & fragment_long == frag]
        if (nrow(temp) > 1){
            temp_reg = lm(apd ~ observed_time_since_infection, data = temp)
            r2 = summary(temp_reg)$r.squared
            slope = coef(temp_reg)['observed_time_since_infection']
            result = data.table(subject_id = sub, fragment = frag, r2 = r2, count = nrow(temp), slope = slope)
            r = rbind(r, result, fill = TRUE)
        }
    }
}

med = r[count > 2, median(r2, na.rm = TRUE), by = fragment]

plot4 = ggplot(r[count > 2]) +
    geom_histogram(aes(x = r2))+
    facet_grid(rows = vars(fragment)) +
    geom_vline(data = med, aes(xintercept = V1), color = 'black', size = 3) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 2)+
    ylab('Observation count')+
    xlab('Coefficient of determination')

name4 = paste0('plots/manuscript_figs/apd_time_r2.pdf')
ggsave(name4, plot = plot4, width = 20, height = 16, units = 'in', dpi = 750, device = cairo_pdf)

vlr = data.table()
tog[fragment == 1, fragment_long := 'Gene region 1\n(within gag)']
tog[fragment != 1, fragment_long := paste0('Gene region ', fragment, '\n(within pol)')]

for (sub in unique(tog$subject_id)){
    for (frag in unique(tog$fragment_long)){
        temp = tog[subject_id == sub & fragment_long == frag & !is.na(vload)]
        if (nrow(temp) > 1){
            temp_reg = lm(apd ~ log10(vload), data = temp)
            r2 = summary(temp_reg)$r.squared
            slope = coef(temp_reg)['log10(vload)']
            result = data.table(subject_id = sub, fragment = frag, r2 = r2, count = nrow(temp), slope = slope)
            vlr = rbind(vlr, result, fill = TRUE)
        }
    }
}

vlmed = vlr[count > 2, median(r2, na.rm = TRUE), by = fragment]

plot4 = ggplot(vlr[count > 2]) +
    geom_histogram(aes(x = r2))+
    facet_grid(rows = vars(fragment)) +
    geom_vline(data = vlmed, aes(xintercept = V1), color = 'black', size = 3) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 30), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 33), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 2)+
    ylab('Observation count')+
    xlab('Coefficient of determination')

name4 = paste0('plots/manuscript_figs/apd_vload_r2.pdf')
ggsave(name4, plot = plot4, width = 20, height = 16, units = 'in', dpi = 750, device = cairo_pdf)


