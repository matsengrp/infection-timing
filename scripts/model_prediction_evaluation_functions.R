configure_newdata <- function(data_path, time_known = FALSE){
    data = check_data(data_path, time_known)
    data = data[!(is.na(pass2_APD))]

    temp_cols = c('ptnum', 'pass2_APD', 'Fragment')
    important_cols = temp_cols[!(temp_cols %in% c('replicate', 'pass2_APD'))]

    subset = data[, ..temp_cols]

    # average APD across replicates
    average_subset = subset[, mean(pass2_APD), by = important_cols]
    setnames(average_subset, 'V1', 'average_APD')

    # convert fragment number to numeric
    average_subset[, fragment_int := as.numeric(substring(Fragment, 2))]
    # assign index to each subject
    average_subset = index_subjects(average_subset)

    data_list = list(
                     observation_count = nrow(average_subset), 
                     subject_count = length(unique(average_subset$ptnum)), 
                     fragment_count = length(unique(average_subset$Fragment)), 
                     apd = average_subset$average_APD, 
                     subject = average_subset$subject_index,
                     subject_id = average_subset$ptnum,
                     fragment = average_subset$fragment_int
                     )
    return(data_list)
}


# simulate y based on new x
predict_posterior <- function(data, model, newdata = TRUE){
    posterior = rstan::extract(model)
    # newdata = ifelse(TRAINING_INFANT_DATA_PATH == TESTING_INFANT_DATA_PATH, FALSE, TRUE)

    if (isTRUE(newdata)){
        stopifnot(all(c('fragment', 'subject_id', 'apd') %in% names(data)))
        predicted_time_since_infection_posteriors = foreach(index = seq(length(data$apd)), .combine = 'cbind') %do% {
            frag = data$fragment[index]
            apd = data$apd[index]
            subject = data$subject_id[index]
            observed_time = data$observed_time[index]
            infection_status = data$infection_status[index]
            if (TIME_CORRECTION_TYPE == 'beta'){
                post_temp = data.table(with(posterior, (fragment_slope_delta_reparameterized[,frag] + baseline_slope + rnorm(subject_slope_delta_mean_estimate, subject_slope_delta_variance_estimate)) * apd))

            } else {
                post_temp = data.table(with(posterior, (fragment_slope_delta[,frag] + baseline_slope + rnorm(subject_slope_delta_mean_estimate, subject_slope_delta_variance_estimate)) * apd))
            }
            colnames(post_temp) = paste0(subject, '_', frag, '_', observed_time, '_', apd)
            post_temp
        }
    } else {
        predicted_time_since_infection_posteriors = posterior$predicted_time_since_infection    
        predicted_time_since_infection_posteriors = data.table(predicted_time_since_infection_posteriors)
        colnames(predicted_time_since_infection_posteriors) = paste0(data$subject_id, '_',data$fragment, '_', data$observed_time, '_', data$apd)
    }

    return(predicted_time_since_infection_posteriors)
}

get_mean_observed_time_to_time_since_infection_prediction <- function(data, model){
    posterior = rstan::extract(model)
    time_correction_posteriors = data.table(posterior$observed_time_to_time_since_infection_correction)
    colnames(time_correction_posteriors) = unique(paste0('subject_', data$subject_id))

    means = time_correction_posteriors[, lapply(.SD, mean)]
    means = means %>% 
        pivot_longer(everything(), names_to = 'subject_id', values_to = 'observed_time_correction') %>%
        as.data.table()

    means$subject_id = substring(means$subject_id, 9)
    return(means)
}

get_mean_prediction <- function(posterior_predictions){
    means = posterior_predictions[, lapply(.SD, mean)]
    means = means %>%
        pivot_longer(everything(), names_to = 'subject', values_to = 'mean_predicted_time_since_infection') %>%
        separate(subject, c('subject_id', 'fragment', 'observed_time', 'apd'), '_') %>%
        as.data.table()
    means[, posterior_variance := lapply(posterior_predictions, var)]
    return(means)
}

predict <- function(data, model){
    newdata = ifelse(TRAINING_INFANT_DATA_PATH == TESTING_INFANT_DATA_PATH, FALSE, TRUE)

    posteriors = predict_posterior(data, model)
    posterior_mean_dt = get_mean_prediction(posteriors)
    #TODO alter this once using complete dataset
    if (isFALSE(newdata)){
        observed_time_corrections = get_mean_observed_time_to_time_since_infection_prediction(data, model)
        posterior_mean_dt = merge(posterior_mean_dt, observed_time_corrections, by = 'subject_id')
        posterior_mean_dt$corrected_observed_time = as.numeric(posterior_mean_dt$observed_time) + as.numeric(posterior_mean_dt$observed_time_correction)
    }

    return(posterior_mean_dt)
}

simulate_apd_time_stan <- function(model){
    apd = seq(0.0001, 0.025, by = 0.0004) 
    post = rstan::extract(model)

    simulations = foreach(index = seq(30), .combine = 'rbind') %do% {
        simulate_subject = rnorm(1,post$subject_slope_delta_mean_estimate[index], post$subject_slope_delta_variance_estimate[index])
        simulate_fragment = rnorm(1,post$fragment_slope_delta_mean_estimate[index], post$fragment_slope_delta_variance_estimate[index])
        baseline_slope = post$baseline_slope[index]
        time_variance = post$time_since_infection_variance_estimate[index]
        time_since_infection = rnorm(length(apd), (simulate_subject + simulate_fragment + baseline_slope) * apd, time_variance)
        data.table(apd = apd, predicted_time_since_infection = time_since_infection)
    }
    return(simulations)
}

get_model_output_filename <- function(input_data){
    path = file.path(OUTPUT_PATH, 'model_predictions')
    dir.create(path, recursive = TRUE)
    data_name = str_remove(input_data, '.csv')
    data_name = str_remove(data_name, '.tsv')
    data_name = str_split(data_name, '/')[[1]][length(str_split(data_name, '/')[[1]])]
    model_name = str_split(MODEL_FILE, '/')[[1]][length(str_split(MODEL_FILE, '/')[[1]])]
    model_name = str_remove(model_name, '.stan')
    file_name = paste0(data_name, '_', model_name, '_results.tsv')
    together = file.path(path, file_name)
    return(together)
}

save_prediction_results <- function(results, input_data){
    file_name = get_model_output_filename(input_data)
    fwrite(results, file_name, sep = '\t')
}


