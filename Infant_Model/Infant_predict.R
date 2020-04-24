###
#' Created on Tue April 20 2020
#'@author: magdalenarussell
###

library(data.table)
library(ggplot2)
library(caret)
library(tidyverse)

source("Infant_functions.R")

#' Data cleaning, APD averaging, Position averaging for testing data
#' 
#' @param data: Dataframe containing testing data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`
#' 
#' @return data_avg: Data.table containing the following columns: "sample", "time", "fragment", "avg_apd1", "avg_apd5", "avg_apd10", "pos_start", and "pos_end"
#' 
data_clean_new_data <- function(data, fragment){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10), mean(HXB2nt_start), mean(HXB2nt_end)), by = .(Sample, ActualTOI..year., Fragment)]
    colnames(data_avg) = c("sample", "time", "fragment", "avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end")
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
predict_infection_time <- function(data, apd_cutoff, fragment, output_filename, type){
    if (type == 'LM_GEE'){
        model = get(paste0(type,"_model", apd_cutoff))
    } else{
        model = get(paste0(type,"_model", apd_cutoff, fragment))
    }
    new_data = data_clean_new_data(data, fragment)
    new_data$estimated_ti = predict(model, newdata = new_data, na.action = na.pass)
    write.csv(new_data, output_filename, row.names = FALSE)
    return(new_data)
}

#' Retrieve regression model object given apd cutoff, fragment, and regression model type
#' 
#' @param apd_cutoff: apd cutoff (options: 1, 5, 10)
#' @param fragment: sequencing region (options: 'F1', 'F2', 'F3')
#' @param type: regression model type (options: "LM", "LM_origin", "LAD", "LAD_origin")
#' 
#' @return model: regression model object
#' 
model_stats <- function(apd_cutoff, fragment, type){
    model = get(paste0(type,"_model", apd_cutoff, fragment))
    return(model)
}