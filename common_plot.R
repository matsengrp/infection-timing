library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)

knitr::opts_chunk$set(fig.width=16, fig.height=8)

#' Organize data by splitting data file by Fragment and by Run number then create new dataframes representing each run for each fragment.
#' 
#' @param data: dataframe containing Fragment and Run columns.
#' @param fragment: column name indicating sequence region names (column entry options are F1, F2, F3).
#' @param run: column name indicating run number for sequencing (column entry options are `Run 1`, `Run 2`, `Run 3`).
#' @return list of dataframes representing each combination of \code{fragment} and \code{run} variables.

make_dataframes <- function(data, fragment, run) {
    data_frags = split(data, fragment)
    F1 <- data_frags$F1
    F2 <- data_frags$F2
    F3 <- data_frags$F3
    data_runs = split(data, run)
    Run1 <- data_runs$`Run 1`
    Run2 <- data_runs$`Run 2`
    Run3 <- data_runs$`Run 3`
    data_split = split(data, list(fragment, run))
    Run1F1 <- data_split$`F1.Run 1`
    Run1F2 <- data_split$`F2.Run 1`
    Run1F3 <- data_split$`F3.Run 1`
    Run2F1 <- data_split$`F1.Run 2`
    Run2F2 <- data_split$`F2.Run 2`
    Run2F3 <- data_split$`F3.Run 2`
    Run3F1 <- data_split$`F1.Run 3`
    Run3F2 <- data_split$`F2.Run 3`
    Run3F3 <- data_split$`F3.Run 3`
    return(list(F1, F2, F3, Run1, Run2, Run3, Run1F1, Run1F2, Run1F3, Run2F1, Run2F2, Run2F3, Run3F1, Run3F2, Run3F3))
}

#' Creates plot title given input variables.
#' 
#' @param apd: integer -- average pairwise diversity value for sequences filtered various ways (options 1, 5, 10)
#' @param fragment: string -- fragment to be plotted (options F1, F2, F3, all)
#' @param run: string -- sequence run to be plotted (options 1, 2, 3, all)
#' @param facet: string -- what variable for the plot to be facetted by (options  Fragment or Run)
#' @param x: string -- what variable to plot on the x axis (options apd or ti or mae)
#' @param y: string -- what variable to plot on the y axis (options apd or eti or vl)
#' @return string to be used as plot title

title <- function(apd, fragment, run, facet1, facet2, x, y) {
    if (y == "eti" && x == "ti"){
        title_name = paste0("Average ETI vs. Actual TI by ", facet1, " for APD-",apd, ", Fragment-", fragment, ", and Run-", run)
    } else if (y == "apd" && x == "ti"){
        title_name = paste0("Average APD vs. Actual TI by ", facet1, " for APD-",apd, ", Fragment-", fragment, ", and Run-", run)
    } else if (x == "vl" && y == "apd"){
        title_name = paste0("Log Viral Load vs. Average APD by Patient for APD-",apd) 
    } else if (x == "vl" && y == "mae"){
        title_name = paste0("Log Viral Load vs. Mean Absolute Error for APD-",apd) 
    } else {
        print("ERROR: incompatible x and y variables")
    }     
    return(title_name)
}

#' Creates plot file name given input variables. 
#' 
#' @param apd: integer -- average pairwise diversity value for sequences filtered various ways (options 1, 5, 10)
#' @param fragment: string -- fragment to be plotted (options F1, F2, F3, all)
#' @param run: string -- sequence run to be plotted (options 1, 2, 3, all)
#' @param facet: string -- what variable for the plot to be facetted by (options  Fragment or Run)
#' @param x: string -- what variable to plot on the x axis (options apd or ti or mae)
#' @param y: string -- what variable to plot on the y axis (options apd or eti or vl)
#' @return string to be used as file name

file <- function(apd, fragment, run, facet1, facet2, x, y) {
    if (y == "eti" && x == "ti"){
        file_name = paste0("AverageETI_ActualTI_by_",facet1,"_for_APD",apd, "_Fragment", fragment, "_Run", run)
    } else if (y == "apd" && x == "ti"){
        file_name = paste0("AverageAPD_ActualTI_by_",facet1,"_for_APD",apd, "_Fragment", fragment, "_Run", run)
    } else if (x == "vl" && y == "apd"){
        file_name = paste0("LogVL_AverageAPD_by_patient","_for_APD",apd) 
    } else if (x == "vl" && y == "mae"){
        file_name = paste0("LogVL_MAE_for_APD",apd) 
    } else {
        print("ERROR: incompatible x and y variables")
    }     
    return(file_name)
}

#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection by fragment for each APD category.
#' 
#' @param data: dataframe containing Fragment and Run columns.
#' @param apd: integer -- average pairwise diversity value for sequences filtered various ways (options 1, 5, 10)
#' @param fragment: string -- fragment to be plotted (options F1, F2, F3, all)
#' @param run: string -- sequence run to be plotted (options 1, 2, 3, all)
#' @param facet1: string -- what variable for the plot to be facetted by (options  Fragment or Run)
#' @param facet2: string -- what variable for the plot to be facetted by (options  Fragment or Run or none)
#' @param points: string -- what variable for the plot points to be shaped by (options Fragment or Run or Sample)
#' @param color_points: string -- what variable for the plot points to be shaped by (options Fragment or Run or Sample)
#' @param x: string -- what variable to plot on the x axis (options ti or apd or mae)
#' @param y: string -- what variable to plot on the y axis (options apd or eti or vl)
#' @return Plot for average ETI vs. actual TI for specified Fragment, Run, APD, and facet

plot_ETI_TI_regression <- function(data, fragment, apd, run, facet1, facet2, points, color_points, x, y) {
    # Set title and file name 
    title_name = title(apd, fragment, run, facet1, facet2, x, y)
    file_name = file(apd, fragment, run, facet1, facet2, x, y)

    # Execute above function to create a list of new dataframes
    dataframe_list = make_dataframes(data, data$Fragment, data$Run)

    # Rename dataframe list entries by Run and Fragment numbers
    names(dataframe_list) = c("F1", "F2", "F3", "Run1", "Run2", "Run3", "Run1F1", "Run1F2", "Run1F3", "Run2F1",  "Run2F2", "Run2F3", "Run3F1", "Run3F2", "Run3F3")

    # Split list of dataframes into individual dataframes
    list2env(dataframe_list ,.GlobalEnv)

    # Assign facetting variables
    if (facet2 == "none"){
        facet_wrap = as.formula(paste("~", facet1))
    } else {
        facet_wrap = as.formula(paste(facet1, "~", facet2))
    }

    # Assign dataframe to use!
    if (fragment == "all" && run == "all"){
        data_arg = data
    } else if (fragment != "all" && run == "all"){
        data_arg = get(fragment)
    } else if (fragment == "all" && run != "all"){
        data_arg = get(paste0("Run",run))
    } else {
        data_arg = get(paste0("Run",run,fragment))
    }

    data_arg$vl_log = log10(data_arg$VL)

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    if (y == "eti" && x == "ti"){
        yvar = paste0("AvgTOI",apd, collapse = NULL)
        ylab = "Average Estimated Time of Infection (years)"

        xvar = "ActualTOI..year."
        xlab = "Actual Time of Infection (years)"

        # Create plot
        p <- ggplot(data=data_arg,aes(data_arg[[yvar]],x=data_arg[[xvar]]))+geom_point(aes(shape=get(paste(points)), color=get(paste(color_points))), size = 2.5) + stat_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + xlim(0,2) + ylim(0,6) + facet_wrap(facet_wrap) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), shape = points, color = color_points)+ stat_regline_equation(label.x = 0, label.y = 5, size = 5) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) 
    } else if (y == "apd" && x == "ti"){
        yvar = paste0("AvgAPD",apd, collapse = NULL)
        ylab = "Average Average Pairwise Diversity (APD)"
    
        xvar = "ActualTOI..year."
        xlab = "Actual Time of Infection (years)"

        # Create plot
        p <- ggplot(data=data_arg,aes(y = data_arg[[yvar]],x=data_arg[[xvar]])) + geom_point(aes(color=get(paste(color_points))), size = 2.5) + facet_wrap(facet_wrap) + geom_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), color = color_points) + stat_cor(method = "pearson", label.x = 1, label.y = 0, size = 5) + theme_bw() + theme(text = element_text(size = 18)) 
    } else if (x == "vl" && y == "apd"){
        xvar = "vl_log"
        xlab = "log10(Viral Load)"

        yvar = paste0("AvgAPD",apd, collapse = NULL)
        ylab = "Average Average Pairwise Diversity (APD)"

        # Create plot
        p <- ggplot(data=data_arg,aes(data_arg[[yvar]],x=data_arg[[xvar]]))+geom_point(aes(color=get(paste(color_points))), size = 2.5) + facet_wrap(facet_wrap) + geom_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), color = color_points) + stat_cor(method = "pearson", label.x = 0, label.y = 0, size = 5)+ theme_bw() + theme(text = element_text(size = 18)) 
    } else if (y == "mae" && x == "vl"){
        xvar = "vl_log"
        xlab = "log10(Viral Load)"

        yvar = paste0("MAE",apd, collapse = NULL)
        ylab = "ETI - TI Mean Absolute Error"

        # Create plot
        p <- ggplot(data=data_arg,aes(data_arg[[yvar]],x=data_arg[[xvar]]))+geom_point(aes(color=get(paste(color_points))), size = 2.5) + facet_wrap(facet_wrap) + geom_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + stat_cor(method = "pearson", label.x = 0, label.y = 0, size = 5) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), color = color_points) + theme_bw() + theme(text = element_text(size = 18)) 
    } else {
        print("ERROR: incompatible x and y variables")
    }

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}

#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection for all fragment and run combinations for each APD status. The shape of the plot points are by patient ID (Sample) or by Run.
#' 
#' @param variable: string -- variable that we are plotting (options apd, eti)
#' @return Lots of plots found in the /plots/ directory!

make_all_plots <- function(x, y){
    # Set working directory to plots folder
    setwd("./plots/")

    # Create a plot for ETI vs. actual TI for APD for all runs facetted by fragment
    if ((x == "vl" && y == "apd") | (x == "vl" && y == "mae")){
        for (diversity in c(1, 5, 10)){
            plot_ETI_TI_regression(all_runs, fragment = "all", apd = diversity, run = "all", facet1 = "Fragment", facet2 = "none",points="none", color_point = "Sample", x = paste(x), y = paste(y))
        }
    } else {
        for (diversity in c(1, 5, 10)){
            fac = "none"
#            for (fac in c("Run", "none")){
                plot_ETI_TI_regression(all_runs, fragment = "all", apd = diversity, run = "all", facet1 = "Fragment", facet2 = fac, points  ="Run",  color_point = "Sample", x = paste(x), y = paste(y))
#            }
        }
    }
}
#' Format all_runs data to contain only patient names, viral load, and actual time of infection
#' 
#' @param variable: data -- dataframe containing columns named `Sample`, `VL`, and `ActualTOI..year.` (among others)
#' @return all_runs_vl -- dataframe containing only columns for `patient`, `time`, and `viral_load`

simplify_df_viral_load <- function(data) {
    if ('Sample' %in% colnames(data)){
        data_simple = data %>% select(Sample, ActualTOI..year., VL)
        colnames(data_simple) = c("patient", "time", "viral_load")
        data_simple_cond = unique(data_simple)
        return(data_simple_cond)
    } else {
        return(data)
    }
}

#' Plot comparison between our data and the Neher test data...specifically patient viral load over time and patient viral load versus apd.
#' 
#' @param variable: string -- variable that we are plotting (options apd, eti)
#' @return Lots of plots found in the /plots/Neher_compare/ directory!

plot_model_compare <- function(data, test_data) {
    setwd("./plots/Neher_compare/")

    data_simple = simplify_df_viral_load(data)
    test_data_simple = simplify_df_viral_load(test_data)
    data_simple$cohort = "Infant"
    test_data_simple$cohort = "Neher"
    data_simple$vl_log = log10(data_simple$viral_load)
    test_data_simple$vl_log = log10(test_data_simple$viral_load)
    both = rbind(data_simple, test_data_simple)
    file_name = "VL_time_Neher_comparison"
    title = "Viral load vs. Time by Sample Cohort"
    xlab = "Time (years)"
    ylab = "log10(Viral Load)"

    p <- ggplot(both, aes(x=time, y=vl_log, color=cohort)) + geom_point(size = 2.5) + geom_smooth(method=lm, se=TRUE, fullrange=TRUE) + ggtitle(title) + labs(x=paste(xlab), y=paste(ylab), color = "Sample Cohort") + theme_bw() + theme(text = element_text(size = 18)) + xlim(0, 2) + stat_cor(method = "pearson", size = 5, label.x = 1.6)

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}