library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)

source("common_regression.R")

knitr::opts_chunk$set(fig.width=16, fig.height=8)

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
        title_name = paste0("Log Viral Load vs. Average Mean Absolute Error for APD-",apd) 
    } else if (x == "none" && y == "none"){
        title_name = paste0("Variable comparison APD",apd) 
    }else if (x == "sp_vl"&& y == "apd_time"){
        title_name = paste0("Set Point Log Viral Load vs. APD Rate of Change for APD-",apd) 
    }else {
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
        file_name = paste0("LogVL_AvgMAE_for_APD",apd) 
    } else if (x == "none" && y == "none"){
        file_name = paste0("Variable_comparison_APD",apd) 
    } else if (x == "sp_vl"&& y == "apd_time"){
        file_name = paste0("SetPtVL_APD_rate_for_APD_",apd) 
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

    # Assign facetting variables
    if (facet2 == "none"){
        facet_wrap = as.formula(paste("~", facet1))
    } else {
        facet_wrap = as.formula(paste(facet1, "~", facet2))
    }

    data_arg = data
    data_arg$vl_log = log10(data_arg$VL)

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    yvar = paste0("AvgTOI",apd, collapse = NULL)
    ylab = "Average Estimated Time of Infection (years)"
    xvar = "ActualTOI..year."
    xlab = "Actual Time of Infection (years)"

    # Create plot
    p <- ggplot(data=data_arg,aes(data_arg[[yvar]],x=data_arg[[xvar]]))+geom_point(aes(shape=get(paste(points)), color=get(paste(color_points))), size = 2.5) + stat_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + xlim(0,2) + ylim(0,6) + facet_wrap(facet_wrap) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), shape = points, color = color_points)+ stat_regline_equation(label.x = 0, label.y = 5, size = 5) + stat_function(fun = func, color = "blue", linetype="dashed", size=1) + theme_bw() + theme(text = element_text(size = 18)) 

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}


#' Plot regressions between different variables (to assess relationships)
#' 
#' @param data: dataframe
#' @param apd: integer -- average pairwise diversity value for sequences filtered various ways (options 1, 5, 10)
#' @param fragment: string -- fragment to be plotted (options F1, F2, F3, all)
#' @param run: string -- sequence run to be plotted (options 1, 2, 3, all)
#' @param facet1: string -- what variable for the plot to be facetted by (options  Fragment or Run)
#' @param facet2: string -- what variable for the plot to be facetted by (options  Fragment or Run or none)
#' @param points: string -- what variable for the plot points to be shaped by (options Fragment or Run or Sample)
#' @param color_points: string -- what variable for the plot points to be shaped by (options Fragment or Run or Sample)
#' @param x: string -- what variable to plot on the x axis (options ti or apd or mae)
#' @param y: string -- what variable to plot on the y axis (options apd or eti or vl)
#' @return Plot for regression bewteen given x and y variables for the given Fragment, Run, APD, and facet

plot_compare_regression <- function(data, fragment, apd, run, facet1, facet2, points, color_points, x, y) {
    # Set title and file name 
    title_name = title(apd, fragment, run, facet1, facet2, x, y)
    file_name = file(apd, fragment, run, facet1, facet2, x, y)

    # Assign facetting variables
    if (facet2 == "none"){
        facet_wrap = as.formula(paste("~", facet1))
    } else {
        facet_wrap = as.formula(paste(facet1, "~", facet2))
    }
    data_arg = data
    data_arg$vl_log = log10(data_arg$VL)

    # Perfect correlation (ETI = Actual TI) function: 
    func <- function(x){x}

    if (y == "apd" && x == "ti"){
        yvar = paste0("AvgAPD",apd, collapse = NULL)
        ylab = "Average Average Pairwise Diversity (APD)"
        ypos1 = 0.021
        ypos2 = 0.022
    
        xvar = "ActualTOI..year."
        xlab = "Actual Time of Infection (years)"
        xpos = 1
    } else if (x == "vl" && y == "apd"){
        xvar = "vl_log"
        xlab = "log10(Viral Load)"
        xpos = 4

        yvar = paste0("AvgAPD",apd, collapse = NULL)
        ylab = "Average Average Pairwise Diversity (APD)"
        ypos1 = 0.04
        ypos2 = 0.038   
    } else if (y == "mae" && x == "vl"){
        xvar = "vl_log"
        xlab = "log10(Viral Load)"
        xpos = 4
        yvar = paste0("AvgMAE",apd, collapse = NULL)
        ylab = "ETI - TI Mean Absolute Error"
        ypos1 = 2.25
        ypos2 = 2.35
    } else if (x == "sp_vl"&& y == "apd_time"){
        data_arg = regress_apd_time(data, apd)
        data_arg$sp_vl_log = log10(data_arg$set_point_VL)
        xvar = "sp_vl_log"
        xlab = "log10(Set Point Viral Load)"
        xpos = 6
        yvar = "apd_time_lm"
        ylab = "APD Rate of Change (diversity/year)"
        ypos1 = 0.017
        ypos2 = 0.018     
    } else {
        print("ERROR: incompatible x and y variables")
    } 
    
    # Create plot
    p <- ggplot(data=data_arg,aes(y = data_arg[[yvar]],x=data_arg[[xvar]])) + geom_point(aes(color=get(paste(color_points))), size = 2.5) + facet_wrap(facet_wrap) + geom_smooth(method="lm", se=TRUE, color= apd, fullrange=TRUE) + ggtitle(title_name) + labs(x=paste(xlab), y=paste(ylab), color =color_points) + stat_cor(method = "pearson", size = 5, label.x = xpos, label.y = ypos1) + stat_regline_equation(label.x = xpos, label.y = ypos2, size = 5) + theme_bw() + theme(text = element_text(size = 18))

    # plot residuals
    #regression = lm(data_arg[[yvar]] ~ data_arg[[xvar]])
    #r <- ggplot(regression, aes(x = .fitted, y = .resid)) + geom_point() + labs(x=paste(xlab), y="Regression Residuals") + ggtitle(paste(title_name, " #Regression Residuals"))

    print(p)
    #print(r)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}

#' Add regression lines by fragment to GGPairs plot...
#' 
#' @param data: dataframe with Fragment column
#' @param fragment: column name indicating sequence region names (column entry options are F1, F2, F3).
#' @return regression lines by Fragment to add to plot
regress_by_frag <- function(data, mapping){
  p <- ggplot(data = data, mapping=mapping) + 
    geom_point() + 
    geom_smooth(method=lm, se=TRUE, fullrange=TRUE, alpha = 0.1)
  p
}

#' Plot scatter matrix to compare relationships between Actual Infection time, APD, and Mean absolute error (MAE)
#' 
#' @param data: dataframe
#' @param apd: integer -- average pairwise diversity value for sequences filtered various ways (options 1, 5, 10)

#' @return Plot for scatter matrix comparing Actual Infection time, APD, and Mean absolute error (MAE)

plot_scatter_matrix <- function(data, apd) {
    data$vl_log = log10(data$VL)
    file_name = file(apd, "none", "none", "none", "none", "none", "none")
    apd_ = paste0("AvgAPD",apd, collapse = NULL)
    mae_ = paste0("AvgMAE",apd, collapse = NULL)
    p <- ggpairs(data=data[,c("ActualTOI..year.", "vl_log", apd_, mae_, "Fragment")], columns = 1:4, title="Actual Infection Time (years) vs. log(Viral Load) and APD", columnLabels = c("Actual TI (years)", "log10(Viral Load)", apd_, mae_), size = 2.5, ggplot2::aes(colour=Fragment, alpha = 0.5),upper = list(continuous = wrap("cor", size=6)), lower = list(continuous = regress_by_frag)) + theme_bw() + theme(text = element_text(size = 18)) 

    print(p)

    path = paste(file_name,".pdf",sep="")
    ggsave(path, plot = last_plot())
}


#' Plot Average Estimated Time of Infection (AvgETI) versus actual time of infection for all fragment and run combinations for each APD status. The shape of the plot points are by patient ID (Sample) or by Run.
#' 
#' @param dataframe: data to plot
#' @param x: string -- variable that we are plotting (options apd, eti)
#' @param y: string -- variable that we are plotting (options apd, eti)
#' @return Lots of plots found in the /plots/ directory!

make_all_plots <- function(data, x, y){
    # Set working directory to plots folder
    setwd("./plots/")

    if (x == "ti" & y == "eti"){
        for (diversity in c(1, 5, 10)){
            plot_ETI_TI_regression(data, fragment = "all", apd = diversity, run = "all", facet1 = "Fragment", facet2 = "none", points  ="Run",  color_point= "Sample", x = paste(x), y = paste(y))
        }
    } else if (x == "none"& y == "none"){
        for (diversity in c(1, 5, 10)){
            plot_scatter_matrix(data, apd = diversity)
        }
    } else {
        for (diversity in c(1, 5, 10)){
            plot_compare_regression(data, fragment = "all", apd = diversity, run = "all", facet1 = "Fragment", facet2 = "none", points = "none",  color_point = "Sample", x = paste(x), y = paste(y))
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
#' @param data -- infant dataframe
#' @param test_data -- dataframe for comparison
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


