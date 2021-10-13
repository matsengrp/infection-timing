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

get_file_path_apd_time_by_subject <- function(){
    path = file.path(PROJECT_PATH, 'plots', 'apd_time', paste0(TIME_CORRECTION_TYPE, '_time_correction'))
    dir.create(path, recursive = TRUE)
    name = paste0(path, '/apd_time_by_subject.pdf')
    return(name)
}

plot_apd_time_by_subject <- function(data){
    data$temp_time = data$mean_predicted_time_since_infection + -1*data$observed_time_correction
    plot = ggplot() +
        geom_point(data = data, aes(x = as.numeric(apd), y = as.numeric(observed_time), color = infection_status), alpha = 0.6, size = 4) +
        geom_smooth(data = data, aes(x = as.numeric(apd), y = temp_time, group = subject_id, color = infection_status), method = 'lm', size = 1, se = FALSE, fullrange = TRUE) +
        facet_grid(~fragment) +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        ylim(-0.25, 2)+
        ylab('Observed (sampling) time') +
        xlab('APD') 
    
    file_name = get_file_path_apd_time_by_subject()
    ggsave(file_name, plot = plot, width = 18, height = 6, units = 'in', dpi = 750, device = cairo_pdf)
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

get_file_path_observed_predicted_time <- function(with_observed_time_correction = FALSE){
    path = file.path(PROJECT_PATH, 'plots', 'observed_predicted_time', paste0(TIME_CORRECTION_TYPE, '_time_correction'))
    dir.create(path, recursive = TRUE)
    if (isFALSE(with_observed_time_correction)){
        name = paste0(path, '/observed_predicted_time.pdf')
    } else {
        name = paste0(path, '/corrected_observed_predicted_time.pdf')
    }
    return(name)
}


plot_observed_predicted_time <- function(data, with_observed_time_correction = FALSE){
    if (isFALSE(with_observed_time_correction)){
        observed_time_temp = 'observed_time'
        xlab = 'Observed Time'
    } else {
        observed_time_temp = 'corrected_observed_time'
        xlab = 'Observed time since infection'
    }    

    plot = ggplot(data) +
        geom_point(aes(x = as.numeric(get(observed_time_temp)), y = as.numeric(mean_predicted_time_since_infection), color = subject_id, shape = fragment), size = 6, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, color = 'black', size = 3) +
        geom_smooth(aes(x = as.numeric(get(observed_time_temp)), y = as.numeric(mean_predicted_time_since_infection)), method = 'lm', size = 3, color = 'red', se = FALSE) +
        theme_cowplot(font_family = 'Arial') +
        theme(axis.text = element_text(size = 20), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 30), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
        background_grid(major = 'xy') +
        xlab(xlab) +
        ylab('Predicted time since infection') +
        labs(color = 'Subject')

    file_name = get_file_path_observed_predicted_time()
    ggsave(file_name, plot = plot, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}
