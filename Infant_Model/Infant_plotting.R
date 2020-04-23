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

#' Creates plot title given input variables.

title <- function(apd, x, y, type, impute) {
    if (y == "eti" && x == "ti"){
        if (impute =='FALSE'){
            title_name = paste0("Average ETI vs. Actual TI for APD-",apd, " using ", type, " Infant Model with impute ", impute) #, '\n', "Actual TI is halfway through third trimester for all")
        } else if (impute == 'TRUE'){
            title_name = paste0("Average ETI vs. Actual TI for APD-",apd, " using ", type, " Infant Model with impute ", impute) #, '\n', "Actual TI is estimated using APD change over time")
        }
    } else {
        print("ERROR: incompatible x and y variables")
    }   
    return(title_name)
}

#' Creates plot file name given input variables. 

file <- function(apd, x, y, type, impute) {
    if (y == "eti" && x == "ti"){
        if (impute =='FALSE'){
            file_name = paste0("AverageETI_ActualTI_for_APD",apd, "_", type, "_Infant_Model_midpoint_actual_TI")
        } else if (impute == 'TRUE'){
            file_name = paste0("AverageETI_ActualTI_for_APD",apd,  "_", type, "_Infant_Model_imputed_actual_TI")
        }
    } else {
        print("ERROR: incompatible x and y variables")
    }     
    return(file_name)
}

#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection by fragment for each APD category.

plot_ETI_TI_regression <- function(train_data, new_data, apd, type, impute) {
    # Set title and file name 
    title_name = title(apd, 'ti', 'eti', type, impute)
    file_name = file(apd,'ti', 'eti', type, impute)

    actual_data = data_clean(train_data, impute)
    # Subset data by indicated fragment
    estimate_together = data.table()
    for (frag in c('F1', 'F2','F3')){
        estimated_data = predict_infection_time(new_data, apd, frag, '../_ignore/train_data_predictions.csv', type, impute)
        estimate_together = rbind(estimate_together, estimated_data)
    }

    together = actual_data[estimate_together, on = c('sample', 'fragment', 'time', 'avg_apd1', 'avg_apd5', 'avg_apd10', 'pos_start', 'pos_end')]    
    # Create plot
    avg_apd = paste0("avg_apd", apd)
    if (type == 'LM'){
        form = 'y ~ x'
    } else if (type == 'LM_log'){
        form = 'y ~ log(x)'
    } else if (type == 'LM_poly'){
        form = 'y ~ poly(x, 3, raw = TRUE)'
    } else if (type == 'LM_origin'){
        form = 'y ~ x + 0'
    } else {
        form = NULL
    }
    #
    #w <- ggplot(data=together,aes(y = actual_ti,x= together[[avg_apd]]))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + #ggtitle(paste0("APD vs. Actual TI for APD", apd, " ", type, " Regression")) + geom_smooth(method = "lm", formula = form, color= apd) + stat_regline_equation#(formula = form, size = 5, label.y = 2.5) + stat_cor(method = "pearson", size = 5, label.y = 2.25) +labs(y='Actual Infection Time', x= #'APD', shape = "Sequence Fragment", color = "Patient")+ theme_bw() + theme(text = element_text(size = 18)) 

    #print(w)
    #path = paste(paste0("ActualTI_vs_APD_for_APD",apd, "_Infant_Model"),".pdf",sep="")
    #ggsave(path, plot = last_plot())

    for (frag in c('F1', 'F2', 'F3')){
        model = model_stats(apd = apd, fragment = frag, type = type, impute = FALSE)
        print(paste0(type, " Model for APD", apd, " and Fragment ", frag, " Summary Statistics:"))
        if (type == 'LM' | type == 'LM_origin'){
            print(model$finalModel)
        } else if (type == 'LAD' | type == 'LAD_origin'){
            print(model)
        }
        print(summary(model))
    }

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    # Create plot
    p <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_regline_equation(size = 5, label.y = 2.5) + stat_cor(method = "pearson", size = 5, label.y = 2.25) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())

    a <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti, group = sample, color = sample))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=FALSE, fullrange=FALSE) + ggtitle(title_name) + labs(x='Actual time since infection', y= 'Estimated time since Infection', shape = "Sequence Fragment", color = "Patient")+ stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) + xlim(-0.25, 2.75) + ylim(-0.25, 2.75)

    print(a)

    together$time_difference = together$estimated_ti - together$actual_ti
    z <- ggplot(together, aes(x=time_difference, fill=fragment)) + facet_wrap(~fragment) + geom_histogram(alpha=0.8, position = 'identity', bins = 10) + labs(x = "ETI-TI [years]", fill = "Fragment") + ggtitle(paste0("Distribution of the Estimation Error by APD", apd, " Model"))+ theme_bw() +  theme(text = element_text(size = 18))
    
    print(z)

    together$time_difference_abs = abs(together$time_difference)
    together2 = together[, mean(time_difference_abs), by = .(actual_ti, fragment)]
    colnames(together2) = c("actual_ti", "fragment", "mean_time_difference_abs")
    x <- ggplot(together2, aes(x=actual_ti, y = mean_time_difference_abs, group=fragment, color = fragment)) + geom_line(size=1) + labs(x = "TI [years]", color = "Fragment", y = "|ETI-TI| [years]") + ggtitle(paste0("Estimation Error Dependence on Time Since Infection (TI) by APD", apd, " Model"))+ theme_bw() +  theme(text = element_text(size = 18))
    
    print(x)
}