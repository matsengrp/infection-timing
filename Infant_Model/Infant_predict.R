###
#' Created on Tue April 20 2020
#'@author: magdalenarussell
###

library(data.table)
library(ggplot2)
library(caret)
library(tidyverse)

source("Infant_functions.R")

# Data cleaning, APD averaging, Position averaging
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
    }
    return(data_avg)
}

## ADD CAPABILITY TO JUST SUBMIT the FOLLWoing!!!!
# needs apd with cutoff of 1, 5, or 10, fragment
predict_infection_time <- function(data, apd_cutoff, fragment, output_filename, type){
    model = get(paste0(type,"_model", apd_cutoff, fragment))
    new_data = data_clean_new_data(data, fragment)
    new_data$estimated_ti = predict(model, newdata = new_data, na.action = na.pass)
    write.csv(new_data, output_filename, row.names = FALSE)
    return(new_data)
}

