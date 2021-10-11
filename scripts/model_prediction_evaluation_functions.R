# simulate y based on new x
predict_posterior <- function(data, model){
    posterior = rstan::extract(model)
    newdata = ifelse(TRAINING_INFANT_DATA_PATH == TESTING_INFANT_DATA_PATH, FALSE, TRUE)

    if (isTRUE(newdata)){
        stopifnot(all(c('fragment', 'subject_id', 'apd') %in% names(data)))
        predicted_time_since_infection_posteriors = foreach(index = seq(length(data$apd)), .combine = 'cbind') %do% {
            frag = data$fragment[index]
            apd = data$apd[index]
            subject = data$subject_id[index]
            observed_time = data$observed_time[index]
            infection_status = data$infection_status[index]
            post_temp = data.table(with(posterior, (fragment_slope_delta[,frag] + baseline_slope) * apd))
            colnames(post_temp) = paste0(infection_status, '_', subject, '_', observed_time, '_', apd, '_', frag)
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
