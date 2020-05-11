library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("multilevel_infant_prediction.R")
source("multilevel_infant_functions.R")

infant_data = read.csv("../_ignore/AllRunsAvg.csv")
setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/models")

models = c("multilevel_subject_var_int_model", "multilevel_subject_fragment_var_int_model", "multilevel_subject_var_int_fragment_var_slope_model", "multilevel_subject_var_int_subject_var_slope_model", "multilevel_subject_var_int_fragment_subject_var_slope_model")
types = c("subject_varying_intercepts", "subject_fragment_varying_intercepts", "subject_varying_intercepts_fragment_varying_slopes", "subject_varying_intercepts_slopes", "subject_varying_intercepts_fragment_subject_varying_slopes")

index = 0
for (model in models){
    index = index + 1
    assign(paste0(model), readRDS(paste0(model,".rds")))
    assign(paste0("posterior_", model), extract.samples(get(model)))
    assign(paste0("results_", model), compute_predictions(indexed_infant_data_cleaned, types[index], get(paste0("posterior_", model))))
}

plot_ETI_TI <- function(data, model){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    plot(data$pred_time, data$time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), xlim = c(-.1, 2.5), ylab = 'Observed Time since infection', xlab = 'Estimated Time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    abline(lm(data$time ~ data$pred_time), col = 'red', lty = 2, lwd = 4)
    legend("topleft", legend=levels(factor(data$subject)), pch=16, col=unique(factor(data$subject)), ncol = 2, cex = 0.75)
}


plot_ETI_TI_by_subject <- function(data, model){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$pred_time, data$time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), xlim = c(-1, 2.5), ylab = 'Observed Time since infection', xlab = 'Estimated Time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    for (samp in unique(data$subject)){
        data_s = data[with(data, subject == samp)]
        abline(lm(data_s$time ~ data_s$pred_time), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
    }
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("topleft", legend=levels(factor(data$subject)), pch=16, col=unique(factor(data$subject)), ncol = 2, cex = 0.75)
}