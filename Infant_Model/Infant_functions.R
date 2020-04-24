###
#' Created on Tue April 20 2020
#' 
#'@author: magdalenarussell
###

library(data.table)
library(ggplot2)
library(caret)
library(tidyverse)
library(L1pack)
library(gee)

source("../common_regression.R")
data = read.csv("../_ignore/AllRunsAvg.csv")


#' Data cleaning, APD averaging, Position averaging for training data
#' 
#' @param data: Dataframe containing training data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' 
#' @return data_avg: Data.table containing the following columns: "sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end", and "actual_ti" (defined as time since infection with infection start as the midpoint of the third trimester)
#' 
data_clean <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10), mean(HXB2nt_start), mean(HXB2nt_end)), by = .(Sample, ActualTOI..year., Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end")
    data_avg$actual_ti = data_avg$time +0.125
    return(data_avg)
}

#' Create lm regression (for indicated fragment and apd) with slope and nonzero y-intercept with apd as regressor using 11-fold cross validation
#' 
#' @param data: Dataframe containing training or testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3')
#' 
#' @return lm regression object for indicated fragment and apd
#' 
create_model_lm <- function(data, cutoff, fragment){
    data_cleaned = data_clean(data)
    
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

#' Create lm regression (for indicated fragment and apd) with slope and zero y-intercept with apd as regressor using 11-fold cross validation
#' 
#' @param data: Dataframe containing training or testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3')
#' 
#' @return lm regression object for indicated fragment and apd with zero y-intercept
#' 
create_model_lm_origin <- function(data, cutoff, fragment){
    data_cleaned = data_clean(data)
    
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

#' Create lad (least absolute deviation) regression (for indicated fragment and apd) with slope and nonzero y-intercept with apd as regressor (here, no 11-fold cross validation)
#' 
#' @param data: Dataframe containing training or testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3')
#' 
#' @return lad regression object for indicated fragment and apd
#' 
create_model_lad <- function(data, cutoff, fragment){
    data_cleaned = data_clean(data)
    
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

#' Create lad (least absolute deviation) regression (for indicated fragment and apd) with slope and zero y-intercept with apd as regressor (here, no 11-fold cross validation)
#' 
#' @param data: Dataframe containing training or testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3')
#' 
#' @return lad regression object for indicated fragment and apd with zero y-intercept
#' 
create_model_lad_origin <- function(data, cutoff, fragment){
    data_cleaned = data_clean(data)
    
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

#' Create lm regression (for indicated apd) with slope and nonzero y-intercept with apd and fragment as regressors using 11-fold cross validation
#' 
#' @param data: Dataframe containing training or testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param cutoff: apd cutoff (options: 1, 5, 10)
#' 
#' @return lm gee regression object for indicated apd
#' 
create_model_lm_gee <- function(data, cutoff, fragment = "all"){
    data_cleaned = data_clean(data)

    # gee requires ordering/sorting by clustering variable
    data_cleaned = data_cleaned[order(sample)]

    # Train the model
    model1 <- gee(actual_ti ~ avg_apd1 + factor(fragment), data = data_cleaned, id = sample, na.action = na.omit)
    model5 <- gee(actual_ti ~ avg_apd5 + factor(fragment), data = data_cleaned, id = sample, na.action = na.omit)
    model10 <- gee(actual_ti ~ avg_apd10 + factor(fragment), data = data_cleaned, id = sample,na.action = na.omit)

    # Summarize the results
    return(get(paste0("model", cutoff)))
}

# Create models for each fragment, apd, model type combination:

for (cut in c(1, 5, 10)){
    assign(paste0('LM_GEE_model', cut), create_model_lm_gee(data, cut))
    for (frag in c('F1', 'F2', 'F3')){
        assign(paste0('LM_model', cut, frag), create_model_lm(data, cut, frag))
        assign(paste0('LM_origin_model', cut, frag), create_model_lm_origin(data, cut, frag))
        assign(paste0('LAD_model', cut, frag), create_model_lad(data, cut, frag))
        assign(paste0('LAD_origin_model', cut, frag), create_model_lad_origin(data, cut, frag))
    }
}


#gee_model_p_val <- function(apd){
#    model = get(paste0('LM_GEE_model', apd))    
#    # Calculate regressor p-values  
#    beta = model$coef # Extract coefficients
#    se = sqrt(diag(model$robust.variance)) # Extract covariance matrix, take diagonal elements, and square root
#    p = 2*(1-pnorm(abs(beta)/se)) 
#    return(p)
#}

