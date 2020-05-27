library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("multilevel_infant_prediction.R")
source("multilevel_infant_functions.R")

infant_data = read.csv("../_ignore/AllRunsAvg.csv")
setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/models")

models = c("multilevel_fragment_subject_var_slope_model_new_no_int", "multilevel_fragment_subject_var_slope_model_time_error", "multilevel_fragment_subject_var_slope_model_time_error2")

types = c("fragment_subject_varying_slopes_no_int", "time_error", "time_error")

models = c("multilevel_fragment_subject_var_slope_model_new2", "multilevel_fragment_subject_var_slope_model_new_no_int", "multilevel_fragment_subject_var_slope_model_time_error2", "multilevel_fragment_subject_var_slope_model_time_error2_var_int")

types = c("fragment_subject_varying_slopes", "fragment_subject_varying_slopes_no_int", "time_error", "time_error_int")

index = 0
for (model in models){
    index = index + 1
    assign(paste0(model), readRDS(paste0(model,".rds")))
    assign(paste0("posterior_", model), extract.samples(get(model)))
    assign(paste0("results_", model), compute_predictions_existing_data(indexed_infant_data_cleaned, types[index], get(paste0("posterior_", model))))
}

plot_APD_TI_no_sim <- function(data, model){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Time (years)', main = paste0('APD versus Time by Patient for Fragment ', model), xlim = c(0,0.03), ylim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$apd, data_s$time, col = col[[samp]], pch=19, cex = 1.25)
            abline(lm(data_s$time ~ data_s$apd), col = col[[samp]], lwd = 2.5)
        }
    legend("topleft", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}

plot_APD_TI_no_sim2 <- function(data, model){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Time (years)', main = paste0('APD versus Time by Patient for Fragment ', model), ylim = c(0,0.03), xlim = c(0, 2.05), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$time, data_s$apd, col = col[[samp]], pch=19, cex = 1.25)
            #abline(lm(data_s$time ~ data_s$apd), col = col[[samp]], lwd = 2.5)
        }
    if (model == "1"){
        abline(a = 0, b = 104.186, lty = 2, lwd = 4)
        text(0.025, 0.1, paste("y = 104.186 x"))
    } else if (model == "2"){
        abline(a = 0, b = 86.027, lty = 2, lwd = 4)
        text(0.025, 0.1, paste("y = 86.027 x"))
    } else if (model == "3"){
        abline(a = 0, b = 130.84, lty = 2, lwd = 4)
        text(0.025, 0.1, paste("y = 130.84 x"))
    } else if (model == "ALL"){
        abline(lm(data$apd ~ 0+data$time), lty = 2, lwd = 4)
    }
    legend("topleft", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}

plot_APD_TI <- function(data, model){
    # the following is great!!!!
    apd = seq(0.0001, 0.03, by = 0.0003)
    simulation = c()
    for (j in apd){
        for (i in 1:50) simulation = c(simulation, simulate_apd_time(model, i, j))
    }
    apd_sim = rep(apd, each = 50)
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(apd_sim, simulation, col = alpha("black", 0.075), pch = 1, xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Actual Time Since Infection (years)', main = paste0('APD versus Actual Time Since Infection by Patient using \n', model), xlim = c(0,0.03), ylim = c(0, 3))
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject ==   samp)]
            points(data_s$apd, data_s$time, col = alpha(col[[samp]], 0.95), pch=19, cex = 1.25)
            abline(lm(data_s$time ~ data_s$apd), col = alpha(col[[samp]], 0.95), lwd = 2.5)
        }
    abline(lm(simulation~apd_sim), col="black", lw = 4)
    legend("topleft", legend=c(levels(factor(data$subject)), "Simulation", "Time ~ APD"), col=c(unique(factor(data$subject)),"black", "black"), lty=c(rep(NA, 12), 1), lwd = c(rep(NA, 12), 2), pch=c(rep(16, 11), 1, NA), ncol = 2, cex = 0.75)
}

plot_APD_TI2 <- function(data, model){
    # the following is great!!!!
    apd = seq(0.0001, 0.03, by = 0.0003)
    simulation = c()
    for (j in apd){
        for (i in 1:50) simulation = c(simulation, simulate_apd_time(model, i, j))
    }
    apd_sim = rep(apd, each = 50)
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(apd_sim, simulation, col = alpha("black", 0.075), pch = 1, xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Actual Time Since Infection (years)', main = paste0('APD versus Actual Time Since Infection by Patient using \n', model), xlim = c(0,0.03), ylim = c(0, 3))
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject ==   samp)]
            points(data_s$apd, data_s$time, col = alpha(col[[samp]], 0.95), pch=19, cex = 1.25)
        }
    abline(lm(simulation~apd_sim), col="black", lw = 4)
    legend("topleft", legend=c(levels(factor(data$subject)), "Simulation", "Time ~ APD"), col=c(unique(factor(data$subject)),"black", "black"), lty=c(rep(NA, 12), 1), lwd = c(rep(NA, 12), 2), pch=c(rep(16, 11), 1, NA), ncol = 2, cex = 0.75)
}

plot_ETI_TI <- function(data, model){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$time, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0,2.5), ylab = 'Estimated time since infection (ETI)', xlab = 'Observed time since infection (OTI)', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    abline(lm(data$pred_time ~ data$time), col = 'red', lty = 2, lwd = 4)
    legend("topleft", legend=c(levels(factor(data$subject)), "ETI = OTI", "ETI ~ OTI"), col=c(unique(factor(data$subject)), "black", "red"), lty=c(rep(NA, 11), 2, 2), lwd = c(rep(NA, 11), 2, 2), pch=c(rep(16, 11), NA, NA), ncol = 2, cex = 0.75)
}


plot_ETI_TI_by_subject <- function(data, model){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$time, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0, 2.5), ylab = 'Estimated time since infection', xlab = 'Observed time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    for (samp in unique(data$subject)){
        data_s = data[with(data, subject == samp)]
        abline(lm(data_s$pred_time ~ data_s$time), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
    }
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("topleft", legend=c(levels(factor(data$subject)), "ETI = OTI"), col=c(unique(factor(data$subject)), "black"), lty=c(rep(NA, 11), 2), lwd = c(rep(NA, 11), 2), pch=c(rep(16, 11), NA), ncol = 2, cex = 0.75)
}