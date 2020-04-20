###
#' Created on Tue April 20 2020
#'@author: magdalenarussell
###

library(data.table)
library(ggplot2)
library(caret)
library(tidyverse)

data = read.csv("../_ignore/AllRunsAvg.csv")

# Data cleaning, APD averaging, Position averaging
data_clean <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10), mean(HXB2nt_start), mean(HXB2nt_end)), by = .(Sample, ActualTOI..year., Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10", "pos_start", "pos_end")
    data_avg$actual_ti = data_avg$time -0.125
    return(data_avg)
}

create_model <- function(data, cutoff, fragment){
    data_cleaned = data_clean(data)
    #cut_arg = paste0('avg_apd', cutoff)
    #cutoffs = c("avg_apd1", "avg_apd5", "avg_apd10")
    #cutoffs2 = cutoffs[!cutoffs == cut_arg]
    #data_cleaned_cutoff = data_cleaned[, !(..cutoffs2)]
    
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
    model5 <- train(actual_ti ~ avg_apd5, data = data_cleaned, method = "lm", trControl = train.control)
    model10 <- train(actual_ti ~ avg_apd10, data = data_cleaned, method = "lm", trControl = train.control)
    
    # Summarize the results
    return(get(paste0("model", cutoff)))
}

model1F1 = create_model(data, 1, 'F1')
model5F1 = create_model(data, 5, 'F1')
model10F1 = create_model(data, 10, 'F1')
model1F2 = create_model(data, 1, 'F2')
model5F2 = create_model(data, 5, 'F2')
model10F2 = create_model(data, 10, 'F2')
model1F3 = create_model(data, 1, 'F3')
model5F3 = create_model(data, 5, 'F3')
model10F3 = create_model(data, 10, 'F3')
