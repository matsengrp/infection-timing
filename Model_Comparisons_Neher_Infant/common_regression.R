library(ggplot2)
library(devtools)
library(gridExtra)
#library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)

knitr::opts_chunk$set(fig.width=16, fig.height=8)


#' Format all_runs data to contain only patient names, viral load, AvgAPD, and actual time of infection
#' 
#' @param data -- dataframe containing columns named `Sample`, `VL`, and `ActualTOI..year.` (among others)
#' @param apd -- diversity cutoff (integer)
#' @param mae -- boolean indicating whether this will be used in an MAE calculation (if so, include AvgTOI)
#' @return dataframe containing only columns for `patient`, `set_point_viral_load`, `fragmetn`,`time`, `actual_viral_load`, and `average_apd`

collapse_df <- function(data, apd, mae) {
    AvgAPD = paste0("AvgAPD",apd, collapse = NULL)
    if (mae == TRUE){
        AvgTOI = paste0("AvgTOI",apd, collapse = NULL)
        data_simple = data %>% select(Sample, Fragment, ActualTOI..year., VL, AvgAPD, AvgTOI, cohort)
        colnames(data_simple) = c("Sample", "Fragment", "ActualTOI..year.", "VL", paste0("AvgAPD", apd), paste0("AvgTOI", apd), "cohort")
    } else if (mae == FALSE){
        data_simple = data %>% select(Sample, Fragment, ActualTOI..year., VL, AvgAPD, cohort)
        colnames(data_simple) = c("Sample", "Fragment", "ActualTOI..year.", "VL", paste0("AvgAPD", apd), "cohort")
    }
    data_simple_cond = unique(data_simple)
    data_simple_cond_dt = data.table(data_simple_cond)
    data_simple_cond_dt_vl = data_simple_cond_dt[VL != "NA",mean(VL)  , by = .(Sample)]
    together = merge(data_simple_cond_dt_vl, data_simple_cond_dt, by = "Sample")
    if (mae == TRUE){
        colnames(together) = c("Sample", "set_point_VL", "Fragment", "ActualTOI..year.", "VL", paste0("AvgAPD", apd), paste0("AvgTOI", apd), "cohort")
    } else if (mae == FALSE){
        colnames(together) = c("Sample", "set_point_VL", "Fragment", "ActualTOI..year.", "VL", paste0("AvgAPD", apd), "cohort")
    }
    return(together)
}

#' Regression function for APD ~ Time by patient, fragment
#' 
#' @param df -- datatable containing columns for `patient`, `fragment`,`time`,  and `average_apd` (produced using `collapse_df()` function)
#' @param frag -- fragment for analysis (options: "F1", "F2", "F3")
#' @param pat -- patient for analysis
#' @param apd -- diversity cutoff (integer)
#' @return linear regression slope and linear regression p-value for relationship between apd and time
regress <- function(df, frag, pat, apd){
    df = data.table(df)
    av_apd = paste0("AvgAPD", apd)
    if (length(df[Sample == pat & Fragment == frag & !is.na(df[[av_apd]])][[av_apd]]) <= 1){
        p = NA
        slope_num = NA
    } else {
        modelobject = lm(df[Fragment == frag & Sample == pat][[av_apd]] ~ df[Fragment == frag & Sample == pat]$ActualTOI..year., na.action = na.exclude)
        slope = coefficients(modelobject)[-1]
        slope_num = as.double(paste(slope[1]))
        f = summary(modelobject)$fstatistic
        p = pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
    }
    return(c(slope_num, p))
}


#' Regression for APD ~ Time by patient, fragment
#' 
#' @param data -- dataframe containing columns named `Sample`, `VL`, and `ActualTOI..year.` (among others)
#' @param apd -- diversity cutoff (integer)
#' @return dataframe containing only columns for `patient`, `set_point_viral_load`, `fragmetn`,`time`, `actual_viral_load`, and `average_apd`

regress_apd_time <- function(data, apd) {
    data = process_data(data)
    data_collapse_df = collapse_df(data, apd, mae=FALSE)
    data_collapse_df = data.table(data_collapse_df)
    av_apd = paste0("AvgAPD", apd)
    lm_df = NULL
    for (frag in c('F1', 'F2', 'F3')){
        Fragment = frag
        for (pat in as.vector(unlist(unique(data_collapse_df$Sample)))){
            Sample = pat
            set_point_VL = data_collapse_df[Sample == pat]$set_point_VL[1]
            cohort = data_collapse_df[Sample == pat]$cohort[1]
            if (all(is.na(data_collapse_df[Sample == pat & Fragment == frag][[av_apd]])) == TRUE){
                apd_time_lm = NA
                apd_time_p_val = NA
            } else {
                apd_time_lm = regress(data_collapse_df, frag, pat, apd)[1]
                apd_time_p_val = regress(data_collapse_df, frag, pat, apd)[2]
            }  
            lm_df = rbind(lm_df, data.frame(Fragment, Sample, set_point_VL, apd_time_lm, apd_time_p_val, cohort))
        }
    }
    return(lm_df) # Note: if dataframe contains an NA/NaN for apd_time_lm, then there was only one value present for the given fragment, patient (no regression possible) and if dataframe contains an NA/NaN for apd_time_p_val, then either no regression was possible, or only two values were present for the given fragment, patient (so no pvalue possible)
}

#' Add average mean absolute error column to dataframe for each apd cutoff level
#' 
#' @param data -- dataframe containing columns named `Sample`, `VL`, and `ActualTOI..year.` (among others)
#' @return dataframe including an average mean absolute error column to dataframe for each apd cutoff level

calc_mae <- function(data) {
    data_mae_all = NULL
    for (apd in c(1, 5, 10)){
        data_collapse_df = collapse_df(data, apd, mae = TRUE)
        data_collapse_df = data.table(data_collapse_df)
        AvgTOI = paste0("AvgTOI",apd, collapse = NULL)
        mae = paste0("AvgMAE",apd, collapse = NULL)
        data_collapse_df$abs_difference = abs(data_collapse_df$ActualTOI..year.-data_collapse_df[[AvgTOI]])
        mae_col  = data_collapse_df[abs_difference != "NA", .(mean(abs_difference)), by = .(Sample, Fragment)]
        data_mae_all = cbind(data_mae_all, mae_col)
    }
    data_mae_all_condensed = data_mae_all[,c(1:3, 6, 9)]
    colnames(data_mae_all_condensed) = c("Sample", "Fragment", "AvgMAE1", "AvgMAE5", "AvgMAE10")
    together = merge(data, data_mae_all_condensed, by = c("Sample", "Fragment"))
    return(together)
}
