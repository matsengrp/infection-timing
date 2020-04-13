library(ggplot2)
library(devtools)
library(gridExtra)
library(ggpubr)
library(data.table)
library(Metrics)
library(GGally)

knitr::opts_chunk$set(fig.width=16, fig.height=8)

#' Format all_runs data to contain only patient names, viral load, AvgAPD, and actual time of infection
#' 
#' @param data -- dataframe containing columns named `Sample`, `VL`, and `ActualTOI..year.` (among others)
#' @param apd -- diversity cutoff (integer)
#' @return dataframe containing only columns for `patient`, `set_point_viral_load`, `fragmetn`,`time`, `actual_viral_load`, and `average_apd`

collapse_df <- function(data, apd) {
    AvgAPD = paste0("AvgAPD",apd, collapse = NULL)
    data_simple = data %>% select(Sample, Fragment, ActualTOI..year., VL, AvgAPD)
    colnames(data_simple) = c("patient", "fragment", "time", "viral_load", paste0("average_apd_", apd))
    data_simple_cond = unique(data_simple)
    data_simple_cond_dt = data.table(data_simple_cond)
    data_simple_cond_dt_vl = data_simple_cond_dt[viral_load != "NA",mean(viral_load)  , by = .(patient)]
    together = merge(data_simple_cond_dt_vl, data_simple_cond_dt, by = "patient")
    colnames(together) = c("patient", "set_point_vl", "fragment", "time", "actual_vl", paste0("average_apd_", apd))
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
    av_apd = paste0("average_apd_", apd)
    if (length(df[patient == pat & fragment == frag & !is.na(df_1dt[[av_apd]])][[av_apd]]) <= 1){
        p = NA
        slope_num = NA
    } else {
        modelobject = lm(df[fragment == frag & patient == pat][[av_apd]] ~ df[fragment == frag & patient == pat]$time, na.action = na.exclude)
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
    data_collapse_df = collapse_df(data, apd)
    data_collapse_df = data.table(data_collapse_df)
    av_apd = paste0("average_apd_", apd)
    lm_df = NULL
    for (frag in c('F1', 'F2', 'F3')){
        fragment = frag
        for (pat in as.vector(unlist(unique(data_collapse_df$patient)))){
            patient = pat
            if (all(is.na(data_collapse_df[patient == pat & fragment == frag][[av_apd]])) == TRUE){
                apd_time_lm = NA
                apd_time_p_val = NA
            } else {
                apd_time_lm = regress(data_collapse_df, frag, pat, apd)[1]
                apd_time_p_val = regress(data_collapse_df, frag, pat, apd)[2]
            }  
            lm_df = rbind(lm_df, data.frame(fragment, patient, apd_time_lm, apd_time_p_val))
        }
    }
    return(lm_df) # Note: if dataframe contains an NA/NaN for apd_time_lm, then there was only one value present for the given fragment, patient (no regression possible) and if dataframe contains an NA/NaN for apd_time_p_val, then either no regression was possible, or only two values were present for the given fragment, patient (so no pvalue possible)
}

