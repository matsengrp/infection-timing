library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)
library(ggfortify)

source("Infant_functions.R")
source("Infant_predict.R")

#' Create plot title given apd cutoff, x and y variables, and regression type
#' 
#' @param apd: apd cutoff (options: 1, 5, 10)
#' @param x: x variable to be plotted (so far, option is "ti")
#' @param y: y variable to be plotted (so far, option is "eti")
#' @param type: regression model type (options: "LM", "LM_origin", "LAD", "LAD_origin")
#' 
title <- function(apd, x, y, type) {
    if (y == "eti" && x == "ti"){
        title_name = paste0("Average ETI vs. Actual TI for APD-",apd, " using ", type, " Infant Model")
    } else {
        print("ERROR: incompatible x and y variables")
    }   
    return(title_name)
}

#' Create file name given apd cutoff, x and y variables, and regression type
#' 
#' @param apd: apd cutoff (options: 1, 5, 10)
#' @param x: x variable to be plotted (so far, option is "ti")
#' @param y: y variable to be plotted (so far, option is "eti")
#' @param type: regression model type (options: "LM", "LM_origin", "LAD", "LAD_origin")
#' 
file <- function(apd, x, y, type) {
    if (y == "eti" && x == "ti"){
        file_name = paste0("AverageETI_ActualTI_for_APD",apd, "_", type, "_Infant_Model_midpoint_actual_TI")
    } else {
        print("ERROR: incompatible x and y variables")
    }     
    return(file_name)
}

#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection by fragment for each APD category.
#'  
#' @param train_data: Dataframe containing training data with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`, `VL`
#' @param new_data: Dataframe containing testing data (for prediction) with columns `APD_1`, `APD_5`, `APD_10`, `HXB2nt_start`, `HXB2nt_end`, `Sample`, `ActualTOI..year.`, `Fragment`
#' @param apd: apd cutoff (options: 1, 5, 10)
#' @param type: regression model type (options: "LM", "LM_origin", "LAD", "LAD_origin")
#' 

plot_ETI_TI_regression <- function(train_data, new_data, apd, type) {
    # Set title and file name 
    title_name = title(apd, 'ti', 'eti', type)
    file_name = file(apd,'ti', 'eti', type)

    actual_data = data_clean(train_data)

    # Subset data by indicated fragment
    if (type == "LM_GEE"){
        estimated_together = predict_infection_time(new_data, apd, fragment = "all", '../_ignore/train_data_predictions.csv', type)
    } else{
        estimate_together = data.table()
        for (frag in c('F1', 'F2','F3')){
            estimated_data = predict_infection_time(new_data, apd, frag, '../_ignore/train_data_predictions.csv', type)
            estimate_together = rbind(estimate_together, estimated_data)
        }
    }

    together = actual_data[estimate_together, on = c('sample', 'fragment', 'time', 'avg_apd1', 'avg_apd5', 'avg_apd10', 'pos_start', 'pos_end')]    
    # Create plot
    avg_apd = paste0("avg_apd", apd)
    if (type == 'LM'){
        form = 'y ~ x'
    } else if (type == 'LM_origin'){
        form = 'y ~ x + 0'
    } else {
        form = NULL
    }
    
    # Plot regression between apd and time since infection 

    #w <- ggplot(data=together,aes(y = actual_ti,x= together[[avg_apd]]))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + #ggtitle(paste0("APD vs. Actual TI for APD", apd, " ", type, " Regression")) + geom_smooth(method = "lm", formula = form, color= apd) + stat_regline_equation#(formula = form, size = 5, label.y = 2.5) + stat_cor(method = "pearson", size = 5, label.y = 2.25) +labs(y='Actual Infection Time', x= #'APD', shape = "Sequence Fragment", color = "Patient")+ theme_bw() + theme(text = element_text(size = 18)) 

    #print(w)
    #path = paste(paste0("ActualTI_vs_APD_for_APD",apd, "_Infant_Model"),".pdf",sep="")
    #ggsave(path, plot = last_plot())

    # Print model statistics: 

    #for (frag in c('F1', 'F2', 'F3')){
    #    model = model_stats(apd = apd, fragment = frag, type = type)
    #    print(paste0(type, " Model for APD", apd, " and Fragment ", frag, " Summary Statistics:"))
    #    if (type == 'LM' | type == 'LM_origin'){
    #        print(model$finalModel)
    #    } else if (type == 'LAD' | type == 'LAD_origin'){
    #        print(model)
    #    }
    #    print(summary(model))
    #}

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    # Create plot comparing estimated time since infection and actual time since infection with a regression line for all patients together
    p <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_regline_equation(size = 5, label.y = 2.5) + stat_cor(method = "pearson", size = 5, label.y = 2.25) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())


    # Create plot comparing estimated time since infection and actual time since infection with regression lines by patient
    a <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti, group = sample, color = sample))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=FALSE, fullrange=FALSE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(a)


    # Create density plot comparing the Distribution of the Estimation Error by APD for each fragment
    together$time_difference = together$estimated_ti - together$actual_ti
    z <- ggplot(together, aes(x=time_difference, fill=fragment)) + facet_wrap(~fragment) + geom_histogram(alpha=0.8, position = 'identity', bins = 10) + labs(x = "ETI-TI [years]", fill = "Fragment") + ggtitle(paste0("Distribution of the Estimation Error by APD", apd, " Model"))+ theme_bw() +  theme(text = element_text(size = 18))
    
    print(z)

    # Create line plot comparing the Estimation Error Dependence on Time Since Infection (TI) by APD for each fragment
    together$time_difference_abs = abs(together$time_difference)
    together2 = together[, mean(time_difference_abs), by = .(actual_ti, fragment)]
    colnames(together2) = c("actual_ti", "fragment", "mean_time_difference_abs")
    x <- ggplot(together2, aes(x=actual_ti, y = mean_time_difference_abs, group=fragment, color = fragment)) + geom_line(size=1) + labs(x = "TI [years]", color = "Fragment", y = "|ETI-TI| [years]") + ggtitle(paste0("Estimation Error Dependence on Time Since Infection (TI) by APD", apd, " Model"))+ theme_bw() +  theme(text = element_text(size = 18))
    
    print(x)
}
