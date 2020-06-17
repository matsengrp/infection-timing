library(rstan)



## PREDICTIONS with training data
# create own link
posterior_link_existing_data_stan <- function(apd, fragment, subject, posterior){
    time <- with(posterior, (subject_slope_delta[,subject] + fragment_slope_delta[,fragment] + baseline_slope) * apd) # here, time is time since infection
    return(time)
}

compute_predictions_existing_data_stan <- function(model_data, posterior){
    prediction_infant_data_cleaned = data.frame()
    for (samp in unique(model_data$subject_index)){
        infant_data_cleaned1 = NULL
        infant_data_cleaned1 = model_data[with(model_data, subject_index == samp)]
        for (frag in unique(infant_data_cleaned1$frag_index)){
            infant_data_cleaned2 = NULL
            apd = NULL
            infant_data_cleaned2 = infant_data_cleaned1[with(infant_data_cleaned1, frag_index == frag)]
            apd = infant_data_cleaned2$apd
            pred.raw <- sapply(1:length(apd), function(i) posterior_link_existing_data_stan(apd[i], frag, samp, posterior))
            pred.p <- apply(pred.raw, 2, mean)
            pred.p.PI <- apply(pred.raw, 2, PI)
            infant_data_cleaned2$pred_time = pred.p
            prediction_infant_data_cleaned = rbind(prediction_infant_data_cleaned, infant_data_cleaned2)
        }   
    }
    return(prediction_infant_data_cleaned)
}

simulate_apd_time_stan = function(model, i, apd){
    post <- extract(get(model))
    simulate_subject = rnorm(1,post$subject_slope_delta_mean_estimate[i], post$subject_slope_delta_variance_estimate[i])
    simulate_fragment = rnorm(1,post$fragment_slope_delta_mean_estimate[i], post$fragment_slope_delta_variance_estimate[i])
    time_since_infection = (simulate_subject + simulate_fragment + post$baseline_slope[i]) * apd
    return(time_since_infection)
}