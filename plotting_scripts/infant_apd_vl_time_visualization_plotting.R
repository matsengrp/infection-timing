library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("multilevel_infant_functions.R")

plot_APD_time <- function(data, fragment){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    if (fragment == 'F1'){
        data = data[fragment == 'F1', mean(apd), by = .(subject, time, fragment, vload)]
    } else if (fragment == 'F2'){
        data = data[fragment == 'F2', mean(apd), by = .(subject, time, fragment, vload)]
    } else if (fragment == 'F3'){
        data = data[fragment == 'F3', mean(apd), by = .(subject, time, fragment, vload)]
    }
    colnames(data) = c('subject', 'time', 'fragment', 'vload', 'apd')
    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Time (years)', main = paste0('APD versus Time by Patient for Fragment ', fragment), xlim = c(0,0.02), ylim = c(0, 2.02), panel.first = grid())
    for (samp in unique(data$subject)){
        data_s = data[with(data, subject == samp)]
        lines(data_s$apd, data_s$time, col = col[[samp]], pch=19, cex = 1.25, type = "b", lty = 1, lwd = 2)
    }
}


plot_VL_time <- function(data){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    data = data[, mean(vload), by = .(subject, time)]
    colnames(data) = c('subject', 'time', 'vload')
    plot(0, type = 'n', ylab = 'Viral Load', xlab = 'Time (years)', main = paste0('Viral Load versus Time by Patient'), ylim = c(0,5007150), xlim = c(0, 2.02), panel.first = grid())
    for (samp in unique(data$subject)){
        data_s = data[with(data, subject == samp)]
        lines(data_s$time, data_s$vload, col = col[[samp]], pch=19, cex = 1.25, type = "b", lty = 1, lwd = 2)
    }
}

plot_VL_apd <- function(data, fragment){
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    if (fragment == 'F1'){
        data = data[fragment == 'F1', mean(apd), by = .(subject, time, fragment, vload)]
    } else if (fragment == 'F2'){
        data = data[fragment == 'F2', mean(apd), by = .(subject, time, fragment, vload)]
    } else if (fragment == 'F3'){
        data = data[fragment == 'F3', mean(apd), by = .(subject, time, fragment, vload)]
    }
    colnames(data) = c('subject', 'time', 'fragment', 'vload', 'apd')
    plot(0, type = 'n', xlab = 'Viral Load', ylab = 'Average Pairwise Diversity (APD1)', main = paste0('Viral Load versus APD by Patient for Fragment ', fragment), xlim = c(0,5007150), ylim = c(0, 0.03), panel.first = grid())
    for (samp in unique(data$subject)){
        data_s = data[with(data, subject == samp)]
        lines(data_s$vload, data_s$apd, col = col[[samp]], pch=19, cex = 1.25, type = "b", lty = 1, lwd = 2)
    }
}
