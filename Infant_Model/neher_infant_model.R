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
library(devtools)
library(gridExtra)
#library(ggpubr)
library(Metrics)
library(GGally)
library(ggfortify)
library(insight)
library(lme4)
library(ggiraphExtra)

source("Infant_functions.R")
source("Infant_predict.R")

source("../Model_Comparisons_Neher_Infant/common_regression.R")

data = read.csv("../_ignore/AllRunsAvg.csv")
neher_data = read.csv("../_ignore/Neher_data_apd_vl.csv")
colnames(neher_data) = c("Sample", "Time", "Fragment", "APD_1", "APD_5", "APD_10", "VL")

#' Data cleaning, APD averaging, Position averaging for training data
#' 
#' @param data: Dataframe containing training data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' 
#' @return data_avg: Data.table containing the following columns: "sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end", and "actual_ti" (defined as time since infection with infection start as the midpoint of the third trimester)
#' 
data_clean_neher <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10)), by = .(Sample, Time, Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    data_avg$actual_ti = data_avg$time
    data_avg = data_avg[, setp_vload := mean(vload), by = .(sample, time, fragment, avg_apd1, avg_apd5, avg_apd10, actual_ti)]
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
create_model_lm_neher <- function(data, cutoff, fragment){
    data_cleaned = data_clean_neher(data)
    
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
    model1 <- train(actual_ti ~ avg_apd1+setp_vload, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model5 <- train(actual_ti ~ avg_apd5+setp_vload, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    model10 <- train(actual_ti ~ avg_apd10+setp_vload, data = data_cleaned, method = "lm", trControl = train.control, na.action = na.pass)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

for (cut in c(1, 5, 10)){
    for (frag in c('F1', 'F2', 'F3')){
        assign(paste0('LM_neher_model', cut, frag), create_model_lm_neher(neher_data, cut, frag))
    }
}

#' Data cleaning, APD averaging, Position averaging for testing data
#' 
#' @param data: Dataframe containing testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`
#' 
#' @return data_avg: Data.table containing the following columns: "sample", "time", "fragment", "avg_apd1", "avg_apd5", "avg_apd10", "pos_start", and "pos_end"
#' 
data_clean_new_data_neher <- function(data, fragment){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10)), by = .(Sample, Time, Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment",  "vload", "avg_apd1", "avg_apd5", "avg_apd10")
    data_avg = data_avg[, setp_vload := mean(vload), by = .(sample, time, fragment, avg_apd1, avg_apd5, avg_apd10)]
    if (fragment == "F1"){
        data_avg = data_avg[fragment == "F1"]
    } else if (fragment == "F2"){
        data_avg = data_avg[fragment == "F2"]
    } else if (fragment == "F3"){
        data_avg = data_avg[fragment == "F3"]
    } else if (fragment == 'all') {
        data_avg = data_avg
    }
    return(data_avg)
}

## ADD CAPABILITY TO JUST SUBMIT the FOLLWoing!!!!

#' Predict estimated time since infection given apd, fragment, regression model type, and testing data
#' 
#' @param data: Dataframe containing testing data (for prediction) with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`
#' @param apd_cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3', 'all')
#' @param output_filename: file name and/or path to write the new data containing estimated time since infection
#' @param type: regression model type (options: "LM", "LM_origin", "LAD", "LAD_origin", "LM_GEE")
#' 
#' @return new_data: Data.table containing the following columns: "sample", "time", "fragment", "avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end", and "estimated_ti"
#' 
predict_infection_time_neher <- function(data, apd_cutoff, fragment, output_filename, type, cohort){
    if (cohort == 'neher'){
        new_data = data_clean_new_data_neher(data, fragment)
    } else {
        new_data = data_clean_new_data(data, fragment)
        new_data = new_data[, setp_vload := mean(vload), by = .(sample, time, fragment, avg_apd1, avg_apd5, avg_apd10)]
    }
    if (type == 'LM_GEE'){
        model = get(paste0(type,"_model", apd_cutoff))
        new_data$estimated_ti = predict(model, data = new_data)
    } else{
        model = get(paste0(type,"_model", apd_cutoff, fragment))
        new_data$estimated_ti = predict(model, newdata = new_data, na.action = na.pass)
    }
    write.csv(new_data, output_filename, row.names = FALSE)
    return(new_data)
}

compile_data_neher <- function(train_data, new_data, apd, type, cohort){
    if (cohort == 'neher'){
        actual_data = data_clean_neher(train_data)
    } else {
        actual_data = data_clean(train_data)
        actual_data = actual_data[, setp_vload := mean(vload), by = .(sample, time, fragment, avg_apd1, avg_apd5, avg_apd10, actual_ti)]
    }
    # Subset data by indicated fragment
    if (type == "LM_GEE"){
        estimate_together = predict_infection_time_neher(new_data, apd, fragment = "all", '../_ignore/train_data_predictions.csv', type, cohort)
    } else{
        estimate_together = data.table()
        for (frag in c('F1', 'F2','F3')){
            estimated_data = predict_infection_time_neher(new_data, apd, frag, '../_ignore/train_data_predictions.csv', type, cohort)
            estimate_together = rbind(estimate_together, estimated_data)
        }
    }
    together = actual_data[estimate_together, on = c('sample', 'fragment', 'time', 'avg_apd1', 'avg_apd5', 'avg_apd10', 'setp_vload', 'vload')]  
    return(together)
}     