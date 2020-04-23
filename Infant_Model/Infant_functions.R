###
#' Created on Tue April 20 2020
#'@author: magdalenarussell
###

library(data.table)
library(ggplot2)
library(caret)
library(tidyverse)
library(L1pack)

source("../common_regression.R")
data = read.csv("../_ignore/AllRunsAvg.csv")

#' Regression function for APD ~ Time by patient, fragment (from common_regression.R)
regress_F3_1 <- function(df, pat){
    df = data.table(df)
    if (length(df[sample == pat]$avg_apd1) <= 1){
        p = NA
        slope_num = NA
    } else {
        modelobject = lm(df[sample == pat]$avg_apd1 ~ df[sample == pat]$time, na.action = na.exclude)
        slope = coefficients(modelobject)[-1]
        intercept = coefficients(modelobject)[-2]
        slope_num = as.double(paste(slope[1]))
        int_num = as.double(paste(intercept[1]))
        f = summary(modelobject)$fstatistic
        p = pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
    }
    return(c(slope_num, int_num, p))
}

# Data cleaning, APD averaging, Position averaging
data_clean <- function(data, impute){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10), mean(HXB2nt_start), mean(HXB2nt_end)), by = .(Sample, ActualTOI..year., Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end")
    if (impute == 'TRUE'){
        data_avg_3 = data_avg[fragment == 'F3']
        lm_df = NULL
        for (pat in as.vector(unlist(unique(data_avg_3$sample)))){
            sample = pat
            apd_time_slope = regress_F3_1(data_avg_3, pat)[1]
            apd_time_intercept = regress_F3_1(data_avg_3, pat)[2]
            apd_time_p_val = regress_F3_1(data_avg_3, pat)[3] 
            lm_df = rbind(lm_df, data.frame(sample, apd_time_slope, apd_time_intercept, apd_time_p_val))
        }
        both = merge(data_avg_3, lm_df)
        both = both[, .(calc_infect_time = -(apd_time_intercept)/apd_time_slope), by = .(sample)]
        both = unique(both)
        data_avg = merge(data_avg, both, by = "sample")
        data_avg = data_avg[calc_infect_time < -0.75, calc_infect_time := -0.125]
        data_avg = data_avg[calc_infect_time > 1, calc_infect_time := -0.125]
        data_avg$actual_ti = data_avg$time - data_avg$calc_infect_time
    } else {
        data_avg$actual_ti = data_avg$time +0.125
    }
    return(data_avg)
}

create_model_lm <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Define training control (11-fold cross validation, repeated 3 times)
    folds <- groupKFold(data_cleaned$sample, k = length(unique(data_cleaned$sample))) 
    train.control <- trainControl(method = "cv", index = folds, savePredictions = "all", returnResamp = "all")

    # Train the model
    model1 <- train(actual_ti ~ avg_apd1, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model5 <- train(actual_ti ~ avg_apd5, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model10 <- train(actual_ti ~ avg_apd10, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

create_model_lm_origin <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Define training control (11-fold cross validation, repeated 3 times)
    folds <- groupKFold(data_cleaned$sample, k = length(unique(data_cleaned$sample))) 
    train.control <- trainControl(method = "cv", index = folds, savePredictions = "all", returnResamp = "all")

    # Train the model
    model1 <- train(actual_ti ~ avg_apd1, data = data_cleaned, method = "lm", trControl = train.control, tuneGrid  = expand.grid(intercept = FALSE), na.action = na.pass)
    model5 <- train(actual_ti ~ avg_apd5, data = data_cleaned, method = "lm", trControl = train.control, tuneGrid  = expand.grid(intercept = FALSE), na.action = na.pass)
    model10 <- train(actual_ti ~ avg_apd10, data = data_cleaned, method = "lm", trControl = train.control, tuneGrid  = expand.grid(intercept = FALSE), na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

create_model_lad <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Train the model
    model1 <- lad(actual_ti ~ avg_apd1, data = data_cleaned, na.action = na.pass)
    model5 <- lad(actual_ti ~ avg_apd5, data = data_cleaned, na.action = na.pass)
    model10 <- lad(actual_ti ~ avg_apd10, data = data_cleaned, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

create_model_lad_origin <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Train the model
    model1 <- lad(actual_ti ~ 0 + avg_apd1, data = data_cleaned, na.action = na.pass)
    model5 <- lad(actual_ti ~ 0 + avg_apd5, data = data_cleaned, na.action = na.pass)
    model10 <- lad(actual_ti ~ 0 + avg_apd10, data = data_cleaned, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

create_model_lm_log <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Define training control (11-fold cross validation, repeated 3 times)
    folds <- groupKFold(data_cleaned$sample, k = length(unique(data_cleaned$sample))) 
    train.control <- trainControl(method = "cv", index = folds, savePredictions = "all", returnResamp = "all")
    data_cleaned = data_cleaned[avg_apd1 == 0, avg_apd1:= NA]
    data_cleaned = data_cleaned[avg_apd5 == 0, avg_apd5:= NA]
    data_cleaned = data_cleaned[avg_apd10 == 0, avg_apd10:= NA]
    # Train the model
    model1 <- train(actual_ti ~ log(avg_apd1), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model5 <- train(actual_ti ~ log(avg_apd5), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model10 <- train(actual_ti ~ log(avg_apd10), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

create_model_lm_poly <- function(data, cutoff, fragment, impute){
    data_cleaned = data_clean(data, impute)
    
    # Subset data by indicated fragment
    if (fragment == "F1"){
        data_cleaned = data_cleaned[fragment == "F1"]
    } else if (fragment == "F2"){
        data_cleaned = data_cleaned[fragment == "F2"]
    } else if (fragment == "F3"){
        data_cleaned = data_cleaned[fragment == "F3"]
    }

    # Define training control (11-fold cross validation, repeated 3 times)
    folds <- groupKFold(data_cleaned$sample, k = length(unique(data_cleaned$sample))) 
    train.control <- trainControl(method = "cv", index = folds, savePredictions = "all", returnResamp = "all")

    # Train the model
    model1 <- train(actual_ti ~ poly(avg_apd1, 3, raw = TRUE), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model5 <- train(actual_ti ~ poly(avg_apd5, 3, raw = TRUE), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model10 <- train(actual_ti ~ poly(avg_apd10, 3, raw = TRUE), data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

for (frag in c('F1', 'F2', 'F3')){
    for (cut in c(1, 5, 10)){
        for (imp in c('FALSE', 'TRUE')){
            assign(paste0('LM_model', cut, frag, '_', imp), create_model_lm(data, cut, frag, imp))
            assign(paste0('LM_origin_model', cut, frag, '_', imp), create_model_lm_origin(data, cut, frag, imp))
            assign(paste0('LAD_model', cut, frag, '_', imp), create_model_lad(data, cut, frag, imp))
            assign(paste0('LAD_origin_model', cut, frag, '_', imp), create_model_lad_origin(data, cut, frag, imp))
            #assign(paste0('LM_log_model', cut, frag, '_', imp), create_model_lm_log(data, cut, frag, imp))
            #assign(paste0('LM_poly_model', cut, frag, '_', imp), create_model_lm_poly(data, cut, frag, imp))
        }
    }
}

