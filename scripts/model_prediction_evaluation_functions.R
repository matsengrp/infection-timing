# simulate y based on new x
predict_posterior <- function(data, model){
    posterior = extract(model)
    newdata = ifelse(TRAINING_INFANT_DATA_PATH == TESTING_INFANT_DATA_PATH, TRUE, FALSE)

    if (isTRUE(newdata)){
    #TODO how to handle unseen subjects..??
        # use baseline slope only?? -- currently doing this!
        # draw a random subject
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
        colnames(predicted_time_since_infection_posteriors) = paste0('subject_', data$subject_id, '_', data$observed_time, '_', data$apd)
    }

    return(predicted_time_since_infection_posteriors)
}

get_mean_prediction <- function(posterior_predictions){
    means = posterior_predictions[, lapply(.SD, mean)]
    means = means %>%
        pivot_longer(everything(), names_to = 'subject', values_to = 'mean_predicted_time_since_infection') %>%
        separate(subject, c('infection_status', 'ptnum', 'observed_time', 'apd', 'fragment'), '_') %>%
        as.data.table()
    means[, posterior_variance := lapply(posterior_predictions, var)]
    return(means)
}

predict <- function(data, model, newdata = TRUE){
    posteriors = predict_posterior(data, model, newdata)
    posterior_mean_dt = get_mean_prediction(posteriors)
    return(posterior_mean_dt)
}

simulate_apd_time_stan <- function(model, i, apd){
    post <- extract(get(model))
    simulate_subject = rnorm(1,post$subject_slope_delta_mean_estimate[i], post$subject_slope_delta_variance_estimate[i])
    simulate_fragment = rnorm(1,post$fragment_slope_delta_mean_estimate[i], post$fragment_slope_delta_variance_estimate[i])
    time_since_infection = (simulate_subject + simulate_fragment + post$baseline_slope[i]) * apd
    return(time_since_infection)
}
