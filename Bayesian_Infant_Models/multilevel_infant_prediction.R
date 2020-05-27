library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)


## PREDICTIONS with training data
# create own link
posterior_link_existing_data <- function(apd, frag_index, subject_index, posterior, type){
    if (type == "fragment_subject_varying_slopes"){
        time <- with(posterior, total_intercept + (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd)
    }
    else if (type == "fragment_subject_varying_slopes_no_int"){
        time <- with(posterior, (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd) # here, time is age
    } else if (type == "time_error"){
        time <- with(posterior, (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd) # here, time is time since infection
    } else if (type == "time_error_int"){
        time <- with(posterior, total_intercept + (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd)
    }
    return(time)
}

compute_predictions_existing_data <- function(model_data, type, posterior){
    prediction_infant_data_cleaned = data.frame()
    for (samp in unique(model_data$subject_index)){
        infant_data_cleaned1 = NULL
        infant_data_cleaned1 = model_data[with(model_data, subject_index == samp)]
        for (frag in unique(infant_data_cleaned1$frag_index)){
            infant_data_cleaned2 = NULL
            apd = NULL
            infant_data_cleaned2 = infant_data_cleaned1[with(infant_data_cleaned1, frag_index == frag)]
            apd = infant_data_cleaned2$apd
            pred.raw <- sapply(1:length(apd), function(i) posterior_link_existing_data(apd[i], frag, samp, posterior, type))
            pred.p <- apply(pred.raw, 2, mean)
            pred.p.PI <- apply(pred.raw, 2, PI)
            infant_data_cleaned2$pred_time = pred.p
            prediction_infant_data_cleaned = rbind(prediction_infant_data_cleaned, infant_data_cleaned2)
        }   
    }
    return(prediction_infant_data_cleaned)
}

simulate_apd_time = function(model, i, apd){
    post <- extract.samples(get(model))
    
    if (model == "multilevel_fragment_subject_var_slope_model_new" | model == "multilevel_fragment_subject_var_slope_model_new2" ) {
        simulate_subject = rnorm(1,post$subject_slope_mean[i], post$subject_slope_sd[i])
        simulate_fragment = rnorm(1,post$fragment_slope_mean[i], post$fragment_slope_sd[i])
        time = post$total_intercept[i] + (simulate_subject + simulate_fragment + post$population_avg_slope[i]) * apd
    } else if (model == "multilevel_fragment_subject_var_slope_model_new_no_int") {
        simulate_subject = rnorm(1,post$subject_slope_mean[i], post$subject_slope_sd[i])
        simulate_fragment = rnorm(1,post$fragment_slope_mean[i], post$fragment_slope_sd[i])
        time = (simulate_subject + simulate_fragment + post$population_avg_slope[i]) * apd
    } else if (model == "multilevel_fragment_subject_var_slope_model_time_error2") {
        simulate_subject = rnorm(1,post$subject_slope_mean[i], post$subject_slope_sd[i])
        simulate_fragment = rnorm(1,post$fragment_slope_mean[i], post$fragment_slope_sd[i])
        time = (simulate_subject + simulate_fragment + post$population_avg_slope[i]) * apd
    } else if (model == "multilevel_fragment_subject_var_slope_model_time_error2_var_int"){
        simulate_subject = rnorm(1,post$subject_slope_mean[i], post$subject_slope_sd[i])
        simulate_fragment = rnorm(1,post$fragment_slope_mean[i], post$fragment_slope_sd[i])
        time = post$total_intercept[i] + (simulate_subject + simulate_fragment + post$population_avg_slope[i]) * apd
    }
    return(time)
}
