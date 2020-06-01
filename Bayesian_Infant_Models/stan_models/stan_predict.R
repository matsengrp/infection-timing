library(rstan)



## PREDICTIONS with training data
# create own link
posterior_link_existing_data_stan <- function(apd, frag_index, subject_index, posterior){
    time <- with(posterior, (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd) # here, time is time since infection
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
    simulate_subject = rnorm(1,post$subject_slope_mean[i], post$subject_slope_sd[i])
    simulate_fragment = rnorm(1,post$fragment_slope_mean[i], post$fragment_slope_sd[i])
    time_since_infection = (simulate_subject + simulate_fragment + post$population_avg_slope[i]) * apd
    return(time_since_infection)
}