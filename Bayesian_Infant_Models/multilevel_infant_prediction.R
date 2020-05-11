library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)


## PREDICTIONS with training data
# create own link
posterior_link <- function(apd, frag_index, subject_index, posterior, type){
    if (type == "subject_varying_intercepts"){
        time <- with(posterior, population_avg_intercept + subject_intercept[subject_index] + population_avg_slope*apd)
    } else if (type == "subject_fragment_varying_intercepts"){
        time <- with(posterior, population_avg_intercept + subject_intercept[subject_index] + fragment_intercept[frag_index] + population_avg_slope*apd)
    } else if (type == "subject_varying_intercepts_fragment_varying_slopes"){
        time <- with(posterior, population_avg_intercept + subject_intercept[,subject_index] + (fragment_slope[,frag_index] + population_avg_slope) * apd)
    } else if (type == "subject_varying_intercepts_slopes"){
        time <- with(posterior, population_avg_intercept + subject_intercept[,subject_index] + (subject_slope[,subject_index] + population_avg_slope) * apd)
    } else if (type == "subject_varying_intercepts_fragment_subject_varying_slopes"){
        time <- with(posterior, population_avg_intercept + subject_intercept[,subject_index] + (subject_slope[,subject_index] + fragment_slope[,frag_index] + population_avg_slope) * apd)
    }
    return(time)
}

compute_predictions <- function(model_data, type, posterior){
    prediction_infant_data_cleaned = data.frame()
    for (samp in unique(model_data$subject_index)){
        infant_data_cleaned1 = NULL
        infant_data_cleaned1 = model_data[with(model_data, subject_index == samp)]
        for (frag in unique(infant_data_cleaned1$frag_index)){
            infant_data_cleaned2 = NULL
            apd = NULL
            infant_data_cleaned2 = infant_data_cleaned1[with(infant_data_cleaned1, frag_index == frag)]
            apd = infant_data_cleaned2$apd
            pred.raw <- sapply(1:length(apd), function(i) posterior_link(apd[i], frag, samp, posterior, type))
            pred.p <- apply(pred.raw, 2, mean)
            pred.p.PI <- apply(pred.raw, 2, PI)
            infant_data_cleaned2$pred_time = pred.p
            prediction_infant_data_cleaned = rbind(prediction_infant_data_cleaned, infant_data_cleaned2)
        }   
    }
    return(prediction_infant_data_cleaned)
}