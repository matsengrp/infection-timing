library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)
library(ggfortify)
library(insight)
library(lme4)
library(ggiraphExtra)

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

plot_ETI_TI_regression <- function(together, apd, type) {
    # Set title and file name 
    title_name = title(apd, 'ti', 'eti', type)
    file_name = file(apd,'ti', 'eti', type)
  
    avg_apd = paste0("avg_apd", apd)


    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    # Create plot comparing estimated time since infection and actual time since infection with a regression line for all patients together
    p <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=TRUE, fullrange=TRUE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_regline_equation(size = 5, label.y = 2.5) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())


    # Create plot comparing estimated time since infection and actual time since infection with regression lines by patient
    a <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti, group = sample, color = sample))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=FALSE, fullrange=FALSE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(a)
}


plot_APD_TI_regression <- function(together, apd, type) {
 
    # Create plot
    avg_apd = paste0("avg_apd", apd)

    if (type == 'LM'){
        form = 'y ~ x'
    } else if (type == 'LM_origin'){
        form = 'y ~ x + 0'
    } else {
        form = NULL
    }
    
    # Print model statistics: 
    if (type == "LM_GEE"){
        model = model_stats(apd = apd, fragment = frag, type = type)
        eqF1 = function(x){coef(model)[2]*x+coef(model)[1]}
        eqF2 = function(x){coef(model)[2]*x+coef(model)[1]+coef(model)[3]}
        eqF3 = function(x){coef(model)[2]*x+coef(model)[1]+coef(model)[4]}
        a = as.numeric(round(coef(model)[1], digits = 3))
        b = as.numeric(round(coef(model)[2], digits = 3))
        c = as.numeric(round(coef(model)[3], digits = 3))
        d = as.numeric(round(coef(model)[4], digits = 3))
        modelF1 = paste0("y = ", a, " + " , b ,"x")
        modelF2 = paste0("y = ", a, " + " , c, " + " ,b ,"x")
        modelF3 = paste0("y = ", a, " + " , d, " + " ,b ,"x")
    } else{
        for (frag in c('F1', 'F2', 'F3')){
            model = model_stats(apd = apd, fragment = frag, type = type)
            if (type == 'LM' | type == 'LM_neher'){
                a = as.numeric(round(coef(model$finalModel)[1], digits = 3))
                b = as.numeric(round(coef(model$finalModel)[2], digits = 3))
            } else if (type == 'LM_origin'){
                a = as.numeric(0)
                b = as.numeric(round(coef(model$finalModel)[1], digits = 3))
            } else if (type == 'LAD'){
                a = as.numeric(round(coef(model)[1], digits = 3))
                b = as.numeric(round(coef(model)[2], digits = 3))
            } else if (type == 'LAD_origin'){
                a = as.numeric(0)
                b = as.numeric(round(coef(model)[1], digits = 3))
            } 
            assign(paste0('model', frag), paste0("y = ", a, " + " , b ,"x"))
        }
    }

    # Plot regression between apd and time since infection (no CI)
    together[fragment == 'F1', labs := paste(modelF1)]
    together[fragment == 'F2', labs := paste(modelF2)]
    together[fragment == 'F3', labs := paste(modelF3)]

    if (type == 'LM'| type == 'LM_origin' | type == 'LM_neher'){
        w <- ggplot(data=together,aes(y = actual_ti,x= together[[avg_apd]]))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + ggtitle(paste0("APD vs. Actual TI for APD", apd, " ", type, " Regression")) + geom_smooth(method = "lm", formula = form)  + geom_text(aes(label = labs), x = 0.015, y = 0.3, size = 5) + labs(y='Actual Time Since Infection', x= paste0('APD', apd), shape = "Fragment", color = "Patient")+ theme_bw() + theme(text = element_text(size = 18)) 
    } else {
        w <- ggplot(data=together,aes(y = actual_ti,x= together[[avg_apd]]))+geom_point(aes(shape=fragment, color = sample), size = 3) + ggtitle(paste0("APD vs. Actual TI for APD", apd, " ", type, " Regression")) + facet_wrap(~fragment) + geom_line(aes(x=together[[avg_apd]], y = estimated_ti), size = 1)+ geom_text(aes(label = labs), x = 0.015, y = 0.3, size = 5)+ labs(y='Actual Time Since Infection', x= paste0('APD', apd) , shape = "Fragment", color = "Patient")+ theme_bw() + theme(text = element_text(size = 18))
    }
    
    print(w)

    path = paste(paste0("ActualTI_vs_APD_for_APD",apd, "_Infant_Model_", type, "_regression"),".pdf",sep="")
    ggsave(path, plot = last_plot())
}

show_model_statistics <- function(apd, type){
    # Print model statistics: 
    if (type == "LM_GEE"){
        model = model_stats(apd = apd, fragment = 'all', type = type)
        print(paste0(type, " Model for APD", apd, " Summary Statistics:"))
        print(summary(model)$coeff)
        p_vals = gee_model_p_val(apd)
        print(p_vals)
    } else{
        for (frag in c('F1', 'F2', 'F3')){
            model = model_stats(apd = apd, fragment = frag, type = type)
            print(paste0(type, " Model for APD", apd, " and Fragment ", frag, " Summary Statistics:"))
            print(summary(model)$coeff)
        }
    }
}


plot_model_evaluation <- function(together, apd, type) {
   
    # Create plot
    avg_apd = paste0("avg_apd", apd)
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