library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)

source("Infant_functions.R")
source("Infant_predict.R")

#' Creates plot title given input variables.

title <- function(apd, x, y) {
    if (y == "eti" && x == "ti"){
        title_name = paste0("Average ETI vs. Actual TI for APD-",apd, " using Infant Model")
    } else {
        print("ERROR: incompatible x and y variables")
    }   
    return(title_name)
}

#' Creates plot file name given input variables. 

file <- function(apd, x, y) {
    if (y == "eti" && x == "ti"){
        file_name = paste0("AverageETI_ActualTI_for_APD",apd, "_Infant_Model")
    } else {
        print("ERROR: incompatible x and y variables")
    }     
    return(file_name)
}

#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection by fragment for each APD category.

plot_ETI_TI_regression <- function(train_data, new_data, apd) {
    # Set title and file name 
    title_name = title(apd, 'ti', 'eti')
    file_name = file(apd,'ti', 'eti')

    actual_data = data_clean(train_data)
    # Subset data by indicated fragment
    estimate_together = data.table()
    for (frag in c('F1', 'F2','F3')){
        estimated_data = predict_infection_time(new_data, apd, frag, '../_ignore/train_data_predictions.csv')
        estimate_together = rbind(estimate_together, estimated_data)
    }
    together = merge(actual_data, estimate_together)

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    # Create plot
    p <- ggplot(data=together,aes(y = estimated_ti,x= actual_ti))+geom_point(aes(shape=fragment, color= sample), size = 2.5) + facet_wrap(~fragment) + stat_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x='Actual Infection Time (Infection starting at third trimester midpoint)', y= 'Estimated Infection Time', shape = "Sequence Fragment", color = "Patient")+ stat_regline_equation(label.x = 1.25, label.y = 2, size = 5) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) 

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}