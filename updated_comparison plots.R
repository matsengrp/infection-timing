library(data.table)
library(dplyr)
library(RColorBrewer)

# Load existing data as `all_runs` file
all_runs = read.csv("_ignore/AllRunsAvg.csv")
all_runs$cohort = "Infant"

# Load Neher paper test data as `Neher_data` file
neher_data = read.csv("_ignore/Neher_data_apd_vl.csv")
neher_data$cohort = "Neher"



# Want to plot time (Age) versus APD, no simulation
plot_APD_time <- function(data, fragment){
    data = as.data.table(data)
    data = data[, subject = .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))
    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'Sampling Time(years)', main = paste0('APD versus Time by Patient for \n Fragment', fragment), xlim = c(0,0.03), ylim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$apd, data_s$time, col = col[[samp]], pch=19, cex = 1.25)
            abline(lm(data_s$time ~ data_s$apd), col = col[[samp]], lwd = 2.5)
    }
    legend("topleft", legend=c(levels(factor(data$Sample))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
}