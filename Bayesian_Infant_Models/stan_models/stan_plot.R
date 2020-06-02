library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("stan_predict.R")

setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models")
source("multilevel_infant_functions.R")

infant_data = read.csv("../_ignore/AllRunsAvg.csv")

setwd("/Users/magdalenarussell/Documents/Matsen_group/infection-timing/Bayesian_Infant_Models/models")


models = c("stan_time_correction_varying_slopes_frag_subject", "stan_time_correction_varying_slopes_frag_subject_more_restricted")

index = 0
for (model in models){
    index = index + 1
    assign(paste0(model), readRDS(paste0(model,".rds")))
    assign(paste0("posterior_", model), extract(get(model)))
    assign(paste0("results_", model), compute_predictions_existing_data_stan(indexed_infant_data_cleaned, get(paste0("posterior_", model))))
}

# Want to plot time (Age) versus APD, no simulation
plot_APD_age_no_sim <- function(data){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Age (years)', main = paste0('APD versus Age by Patient'), xlim = c(0,0.03), ylim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$apd, data_s$time, col = col[[samp]], pch=19, cex = 1.25)
            abline(lm(data_s$time ~ data_s$apd), col = col[[samp]], lwd = 2.5)
        }
    legend("topleft", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}

# Want to plot time_since_infection (predicted by model) versus APD, no simulation
plot_APD_TI_no_sim <- function(data, model){
    assign("data", data)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$apd, data$time_since_infection, xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Time Since Infection (years)', main = paste0('APD versus Time Since Infection by Patient'), xlim = c(0,0.03), ylim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$apd, data_s$time_since_infection, col = col[[samp]], pch=19, cex = 1.25)
        }
    legend("topleft", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}

# Want to plot time_since_infection (predicted by model) versus age
plot_age_TI <- function(data, model){
    assign("data", data)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$time, data$time_since_infection, ylab = 'Time since infection (determined in model fitting)', xlab = "Age (years)", main = paste0('APD versus Time Since Infection by Patient \n Model fit time since infection values \n', model), xlim = c(0, 0.5), ylim = c(0, 0.5), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$time, data_s$time_since_infection, col = col[[samp]], pch=19, cex = 1.25)
            abline(lm(data_s$time_since_infection ~ data_s$time), col = col[[samp]], lwd = 2.5)
        }
    legend("bottomright", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}


# Want to plot time_since_infection (predicted by model) versus age
plot_age_predTI <- function(data, model){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(data$time, data$pred_time, ylab = 'Time Since Infection (predicted using model)', xlab = "Age (years)", main = paste0('APD versus Time Since Infection by Patient \n Model predicted time since infection values \n', model), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$time, data_s$pred_time, col = col[[samp]], pch=19, cex = 1.25)
            abline(lm(data_s$pred_time ~ data_s$time), col = col[[samp]], lwd = 2.5)
        }
    legend("bottomright", legend=c(levels(factor(data$subject))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}


plot_APD_TI <- function(data, model){
    assign("data", data)
    # the following is great!!!!
    apd = seq(0.0001, 0.03, by = 0.0004)
    simulation = c()
    for (j in apd){
        for (i in 1:30) simulation = c(simulation, simulate_apd_time_stan(model, i, j))
    }
    apd_sim = rep(apd, each = 30)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(apd_sim, simulation, pch = 1, col = alpha("black", 0.3), xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Time Since Infection (years) (determined in model fitting)', main = paste0('APD versus Time Since Infection by Patient \n', model), xlim = c(0,0.03), ylim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$apd, data_s$time_since_infection, col = col[[samp]], pch=19, cex = 1.25)
        }
    abline(lm(simulation~apd_sim), col="black", lw = 4)
    cf <- round(coef(lm(simulation~apd_sim)), 2)
    eq <- paste0("time = ", cf[1], "+", cf[2],"* APD")
    mtext(eq, 1, line=-2)
    legend("topleft", legend=c(levels(factor(data$subject)), "Simulation", "Time ~ APD"), col=c(unique(factor(data$subject)),"black", "black"), lty=c(rep(NA, 12), 1), lwd = c(rep(NA, 12), 2), pch=c(rep(16, 11), 1, NA), ncol = 2, cex = 0.75)
}

plot_ETI_TI <- function(data, model, time_type){ 
    assign("data", data)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    if (time_type == "TI"){
        plot(data$time_since_infection, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0,2.5), ylab = 'Estimated time since infection (ETI)', xlab = 'Observed time since infection (OTI)', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
        abline(a = 0, b = 1, lty = 2, lwd = 4)
        abline(lm(data$pred_time ~ data$time_since_infection), col = 'red', lty = 2, lwd = 4)
        cf <- round(coef(lm(data$pred_time ~ data$time_since_infection)), 2)    
    } else if (time_type == "age"){
        plot(data$time, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0,2.5), ylab = 'Estimated time since infection (ETI)', xlab = 'Age at sampling time (years)', main = paste0('Estimated versus Age at sampling time by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
        abline(a = 0, b = 1, lty = 2, lwd = 4)
        abline(lm(data$pred_time ~ data$time), col = 'red', lty = 2, lwd = 4)
        cf <- round(coef(lm(data$pred_time ~ data$time)), 2)    
    }
    eq <- paste0("Estimated Time = ", cf[1], "+", cf[2],"* Observed Time")
    mtext(eq, 1, line=-2)
    legend("topleft", legend=c(levels(factor(data$subject)), "ETI = OTI", "ETI ~ OTI"), col=c(unique(factor(data$subject)), "black", "red"), lty=c(rep(NA, 11), 2, 2), lwd = c(rep(NA, 11), 2, 2), pch=c(rep(16, 11), NA, NA), ncol = 2, cex = 0.75)
}


plot_ETI_TI_by_subject <- function(data, model, time_type){ 
    assign("data", data)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    if (time_type == "TI"){
        plot(data$time_since_infection, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0, 2.5), ylab = 'Estimated time since infection', xlab = 'Observed time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
        for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            abline(lm(data_s$pred_time ~ data_s$time_since_infection), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
        }
    } else if (time_type == "age"){
        plot(data$time, data$pred_time, col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0, 2.5), ylab = 'Estimated time since infection', xlab = 'Age at sampling time (years)', main = paste0('Estimated versus Age at sampling time by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
        for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            abline(lm(data_s$pred_time ~ data_s$time), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
        }
    }
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("topleft", legend=c(levels(factor(data$subject)), "ETI = OTI"), col=c(unique(factor(data$subject)), "black"), lty=c(rep(NA, 11), 2), lwd = c(rep(NA, 11), 2), pch=c(rep(16, 11), NA), ncol = 2, cex = 0.75)
}


plot_ETI_TI_subject_specified <- function(data, model, subject_list){ 
    assign("data", data)
    posterior = get(paste0("posterior_", model))
    times_since_infection = apply(posterior$time_since_infection, 2, mean)
    data$time_since_infection = times_since_infection
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 3, name = 'Set3'))
    col = setNames(palette(), levels(data$fragment))
    plot(0, type = 'n', col = c(factor(data$subject), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2.5), xlim = c(0, 2.5), ylab = 'Estimated time since infection', xlab = 'Observed time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient using \n', model), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    for (samp in subject_list){
        data_s = data[with(data, subject == samp)]
        for (frag in unique(data$frag_index)){
            data_s_f = data_s[with(data_s, frag_index == frag)]
            points(data_s_f$apd, data_s_f$time_since_infection, col = col[[frag]], pch=19, cex = 1.25)
            if (nrow(data_s_f) > 1){
                abline(lm(data_s_f$pred_time ~ data_s_f$time_since_infection), col = col[[frag]], lwd = 2.5)
            }
        }
    }
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    #legend()
}

