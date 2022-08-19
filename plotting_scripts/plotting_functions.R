get_file_path_apd_time <- function(time_type, real_data = NULL){
    path = file.path(PROJECT_PATH, 'plots', 'apd_time', paste0(TIME_CORRECTION_TYPE, '_time_correction'))
    dir.create(path, recursive = TRUE)
    if (is.null(real_data)){
        name = paste0(path, '/apd_', time_type, '.pdf')
    } else {
        name = paste0(path, '/apd_', time_type, '_with_real_data.pdf')
    }
    return(name)
}

get_file_path_apd_observed_time_by_subject <- function(type){
    path = file.path(PROJECT_PATH, 'plots', 'apd_observed_time')
    dir.create(path, recursive = TRUE)
    name = paste0(path, '/apd_time_', type, '.pdf')
    return(name)
}

plot_apd_observed_time_by_subject <- function(data){
    plot = ggplot() +
        geom_point(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time), color = as.factor(subject)), alpha = 0.6, size = 4) +
        geom_smooth(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time), color = as.factor(subject)), method = 'lm', size = 1, se = FALSE, fullrange = TRUE) +
        facet_grid(~fragment) +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        ylim(0, 0.04)+
        xlab('Observed (sampling) time') +
        ylab('APD') +
        guides(fill=guide_legend(title="Participant"))
    
    file_name = get_file_path_apd_observed_time_by_subject('by_subject')
    ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_apd_observed_time_by_subject_color_timing <- function(data){
    plot = ggplot() +
        geom_point(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time), color = as.factor(subject)), alpha = 0.6, size = 4) +
        geom_smooth(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time), linetype = as.factor(infection_status), color = as.factor(subject)), method = 'lm', size = 1, se = FALSE, fullrange = TRUE) +
        facet_grid(~fragment) +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        ylim(0, 0.04)+
        xlab('Observed (sampling) time') +
        ylab('APD') +
        guides(fill=guide_legend(title="Participant"))
    
    file_name = get_file_path_apd_observed_time_by_subject('by_subject_color_timing')
    ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}


get_file_path_apd_time_by_subject <- function(type){
    path = file.path(PROJECT_PATH, 'plots', 'apd_observed_time')
    dir.create(path, recursive = TRUE)
    name = paste0(path, '/apd_time_', type, '.pdf')
    return(name)
}

plot_apd_time_by_subject <- function(data){
    # data$temp_time = data$mean_predicted_time_since_infection + -1*data$observed_time_correction

    nice_fragments = c('Sequencing region 1\n(within gag)', 'Sequencing region 2\n(within pol)', 'Sequencing region 3\n(within pol)')
    data$fragment_long = mapvalues(data$fragment, from = seq(1,3), to = nice_fragments)
    plot = ggplot() +
        # geom_point(data = data, aes(x = as.numeric(apd), y = as.numeric(temp_time), color = infection_status), alpha = 0.6, size = 4) +
        geom_point(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time_since_infection), color = as.factor(subject)), alpha = 0.6, size = 6) +
        # geom_smooth(data = data, aes(x = as.numeric(apd), y = temp_time, group = subject_id, color = infection_status), method = 'lm', size = 1, se = FALSE, fullrange = TRUE) +
        geom_smooth(data = data, aes(y = as.numeric(apd), x = observed_time_since_infection, color = as.factor(subject)), method = 'lm', size = 2.5, se = FALSE, fullrange = TRUE) +
        facet_grid(~fragment_long) +
        theme_cowplot(font_family = 'Arial') +
        theme(legend.position = 'none', axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        ylim(0, 0.04)+
        xlab('Estimated time since infection (years)') +
        ylab('APD') 
    
    file_name = get_file_path_apd_time_by_subject('by_subject')
    ggsave(file_name, plot = plot, width = 18, height = 9, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_apd_time_by_subject_color_timing <- function(data){
    # data$temp_time = data$mean_predicted_time_since_infection + -1*data$observed_time_correction

    nice_fragments = c('Sequencing region 1\n(within gag)', 'Sequencing region 2\n(within pol)', 'Sequencing region 3\n(within pol)')
    data$fragment_long = mapvalues(data$fragment, from = seq(1,3), to = nice_fragments)
    plot = ggplot() +
        # geom_point(data = data, aes(x = as.numeric(apd), y = as.numeric(temp_time), color = infection_status), alpha = 0.6, size = 4) +
        geom_point(data = data, aes(y = as.numeric(apd), x = as.numeric(observed_time_since_infection), color = as.factor(subject)), alpha = 0.6, size = 6) +
        # geom_smooth(data = data, aes(x = as.numeric(apd), y = temp_time, group = subject_id, color = infection_status), method = 'lm', size = 1, se = FALSE, fullrange = TRUE) +
        geom_smooth(data = data, aes(y = as.numeric(apd), x = observed_time_since_infection, linetype = as.factor(infection_status), color = as.factor(subject)), method = 'lm', size = 2.5, se = FALSE, fullrange = TRUE) +
        facet_grid(~fragment_long) +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 25), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 30), axis.line = element_blank(), text = element_text(size = 37), axis.ticks = element_line(color = 'gray60', size = 1.5), legend.key.width = unit(3,"cm")) +
        background_grid(major = 'xy') +
        ylim(0, 0.04)+
        xlab('Sampling time (years)') +
        ylab('APD')+
        labs(linetype = 'Infection time') +
        guides(color = 'none')
    
    file_name = get_file_path_apd_time_by_subject('by_subject_color_timing')
    ggsave(file_name, plot = plot, width = 22, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}


plot_apd_time_all <- function(data, by_subject_regression_sim = NULL){
    plot = ggplot() +
        geom_point(data = data, aes(y = as.numeric(apd), x = observed_time_since_infection, color = as.factor(subject_id)), alpha = 0.6, size = 6) 

    if (!is.null(by_subject_regression_sim)){
        plot = plot + 
            geom_line(data = by_subject_regression_sim, aes(y = predicted_apd, x = observed_time_since_infection, color = as.factor(subject_id)), size = 1.5, alpha = 0.6)
    }

    plot = plot +
        geom_smooth(data = data, aes(y = as.numeric(apd), x = observed_time_since_infection), method = 'lm', size = 3, color = 'black') +
        theme_cowplot(font_family = 'Arial') +
        theme(legend.position = 'none', axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        # scale_x_continuous(breaks=seq(0, 0.03, 0.005)) +
        # scale_y_continuous(breaks=seq(0, 3, 0.5)) +
        ylab('Average pairwise diversity') +
        xlab('Sampling time (years)') 

    file_name = get_file_path_apd_observed_time_by_subject('all')
   
    ggsave(file_name, plot = plot, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}


plot_apd_time <- function(sim_data, data = NULL){
    plot = ggplot() +
        geom_point(data = sim_data, aes(x = as.numeric(apd), y = predicted_time_since_infection), color = 'grey', alpha = 0.3, size = 6) +
        geom_smooth(data = sim_data, aes(x = as.numeric(apd), y = predicted_time_since_infection), method = 'lm', size = 3, color = 'black') +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        # scale_x_continuous(breaks=seq(0, 0.03, 0.005)) +
        # scale_y_continuous(breaks=seq(0, 3, 0.5)) +
        xlab('APD') +
        ylab('Time Since Infection') 
    
    if (!is.null(data)){
        plot = plot +
            geom_point(data = data, aes(x = as.numeric(apd), y = as.numeric(corrected_observed_time), color = as.character(subject_id)), size = 6, alpha = 0.7) +
            labs(color = 'Subject')

    }
    file_name = get_file_path_apd_time('predicted_time_since_infection', data)
    ggsave(file_name, plot = plot, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}

get_file_path_observed_predicted_time <- function(with_observed_time_correction = FALSE, loocv = FALSE, hist = FALSE, adult_data){
    path = file.path(PROJECT_PATH, 'plots', 'observed_predicted_time', paste0(TIME_CORRECTION_TYPE, '_time_correction'))
    dir.create(path, recursive = TRUE)
    if (isFALSE(with_observed_time_correction)){
        name = paste0(path, '/observed_predicted_time')
    } else {
        name = paste0(path, '/corrected_observed_predicted_time')
    }
    if (isFALSE(loocv)){
        name = paste0(name)
    } else {
        name = paste0(name, '_loocv')
    }
    if (isFALSE(adult_data)){
        name = paste0(name)
    } else {
        name = paste0(name, '_ADULT_MODEL')
    }
    if (isTRUE(hist)){
        name = paste0(name, '_hist.pdf')
    } else {
        name = paste0(name, '.pdf')
    }
    return(name)
}


plot_observed_predicted_time <- function(data, with_observed_time_correction = FALSE, loocv = FALSE, with_subject_legend = TRUE, with_fragment_legend = TRUE, adult_data = FALSE, xlimits = NA, ylimits = NA, write_plot = TRUE){
    if (isTRUE(with_observed_time_correction)){
        observed_time_temp = 'corrected_observed_time_since_infection'
        xlab = 'Observed time since infection (corrected)'
    } else {
        observed_time_temp = 'observed_time_since_infection'
        xlab = 'Estimated time since infection'
    }    
    stopifnot(observed_time_temp %in% colnames(data))


    plot = ggplot(data) +
        geom_point(aes(x = as.numeric(get(observed_time_temp)), y = as.numeric(mean_time_since_infection_estimate), color = as.factor(subject_id), shape = as.factor(fragment)), size = 6, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, color = 'black', size = 3) +
        geom_smooth(aes(x = as.numeric(get(observed_time_temp)), y = as.numeric(mean_time_since_infection_estimate)), method = 'lm', size = 3, color = 'red', se = FALSE) +
        theme_cowplot(font_family = 'Avenir') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        xlab(xlab) +
        ylab('Predicted time since infection') +
        labs(color = 'Subject', shape = 'Sequencing\nregion') +
        panel_border(color = 'gray60', size = 2)

    if (isFALSE(with_subject_legend)){
        plot = plot + guides(color = 'none')
    }

    if (isFALSE(with_fragment_legend)){
        plot = plot + guides(shape = 'none')
    }

    if (!(is.na(xlimits))){
        plot = plot + scale_x_continuous(limits = xlimits)
    }

    if (!(is.na(ylimits))){
        plot = plot + scale_y_continuous(limits = ylimits)
    }


    if (isTRUE(write_plot)){
        file_name = get_file_path_observed_predicted_time(with_observed_time_correction, loocv, hist = FALSE, adult_data)
        ggsave(file_name, plot = plot, width = 9, height = 9.5, units = 'in', dpi = 750, device = cairo_pdf)
    }
    return(plot)
}

plot_observed_predicted_time_histogram <- function(data, with_observed_time_correction = FALSE, loocv = FALSE, adult_data = FALSE, xlimits = NA, ylimits = NA, write_plot = TRUE){
    if (isTRUE(with_observed_time_correction)){
        observed_time_temp = 'corrected_observed_time_since_infection'
        xlab = 'Predicted - Observed (years)'
    } else {
        observed_time_temp = 'observed_time_since_infection'
        xlab = 'Predicted - Estimated (years)'
    }    
    stopifnot(observed_time_temp %in% colnames(data))
    
    data$difference = data$mean_time_since_infection_estimate - data[[observed_time_temp]]

    plot = ggplot(data) +
        geom_histogram(aes(x = difference), alpha = 0.7) +
        facet_grid(rows = vars(fragment))+
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        theme_cowplot(font_family = 'Avenir') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        xlab(xlab) +
        ylab('density') +
        panel_border(color = 'gray60', size = 2)
    if (!(is.na(xlimits))){
        plot = plot + scale_x_continuous(limits = xlimits)
    }

    if (!(is.na(ylimits))){
        plot = plot + scale_y_continuous(limits = ylimits)
    }

    if (isTRUE(write_plot)){
        file_name = get_file_path_observed_predicted_time(with_observed_time_correction, loocv, hist = TRUE, adult_data)
        ggsave(file_name, plot = plot, width = 13, height = 11.25, units = 'in', dpi = 750, device = cairo_pdf)
    } 
    return(plot)
}

plot_observed_time_versus_time_diff <- function(data){
    if (isTRUE(with_observed_time_correction)){
        observed_time_temp = 'corrected_observed_time_since_infection'
        xlab = '|Predicted - Observed (years)|'
    } else {
        observed_time_temp = 'observed_time_since_infection'
        xlab = '|Predicted - Estimated (years)|'
    }    
    stopifnot(observed_time_temp %in% colnames(data))
    
    data$difference = data$mean_time_since_infection_estimate - data[[observed_time_temp]]
    plot_data = data[, mean(abs(difference)), by = .(fragment, observed_time_since_infection)]
    plot = ggplot(plot_data) +
        geom_point(aes(x = observed_time_since_infection, y = V1, color = as.factor(fragment))) 
}
