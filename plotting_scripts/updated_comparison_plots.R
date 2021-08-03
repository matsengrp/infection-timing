library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# Load existing data as `all_runs` file
all_runs = read.csv("_ignore/AllRunsAvg.csv")
all_runs$cohort = "Infant"

# Load Neher paper test data as `Neher_data` file
neher_data = read.csv("_ignore/Neher_data_apd_vl.csv")
neher_data$cohort = "Neher"



# Want to plot time (Age) versus APD, no simulation
plot_APD_time <- function(data, fragment, legend){
    data = as.data.table(data)
    data = data[, subject := .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]
    data = data[Fragment == fragment]

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))

    filename = paste0('figures/APD_versus_Time_by_Patient_for_', fragment, '.png')

    png(file=filename, width = 6, height = 5, units = 'in', res = 200)

    plot(0, type = 'n', ylab = 'Average Pairwise Diversity (APD1)', xlab = 'Sampling Time(years)', main = paste0('APD versus Time by Patient for \n Fragment ', fragment), ylim = c(0,0.03), xlim = c(0, 3), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$time, data_s$apd, col = col[[samp]], pch=19, cex = 1.25)
            data_s = data_s[!is.na(time) & !is.na(apd)]

            if (nrow(data_s) >= 3){
                abline(lm(data_s$apd ~ data_s$time), col = col[[samp]], lwd = 2.5)
            }
    }
    if (legend == "True"){
        legend("topleft", legend=c(levels(factor(data$Sample))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
    }
    
    dev.off()
}

plot_APD_time(all_runs, 'F1', legend = 'True')
plot_APD_time(all_runs, 'F2', legend = 'False')
plot_APD_time(all_runs, 'F3', legend = 'False')

# Want to plot time (Age) versus APD, no simulation
plot_APD_vl <- function(data, fragment, legend){
    data = as.data.table(data)
    data = data[, subject := .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]
    data = data[Fragment == fragment]
    data$vl = log10(data$VL)

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))

    filename = paste0('figures/APD_versus_VL_by_Patient_for_', fragment, '.png')

    png(file=filename, width = 6, height = 5, units = 'in', res = 200)

    plot(0, type = 'n', xlab = 'Average Pairwise Diversity (APD1)', ylab = 'log10 Viral Load', main = paste0('APD versus Viral Load for \n Fragment ', fragment), xlim = c(0,0.03), panel.first = grid())
    for (samp in unique(data$subject)){
            data_s = data[with(data, subject == samp)]
            points(data_s$time, data_s$apd, col = col[[samp]], pch=19, cex = 1.25)
    }
    abline(lm(data$apd ~ data$vl), lwd = 2.5)
    if (legend == "True"){
        legend("topleft", legend=c(levels(factor(data$Sample))), col=c(unique(factor(data$subject))), lty=c(rep(NA, 11)), lwd = c(rep(NA, 11)), pch=c(rep(19, 11)), ncol = 2, cex = 0.75)
    }
    
    dev.off()
}

plot_APD_vl(all_runs, 'F1', legend = 'True')
plot_APD_vl(all_runs, 'F2', legend = 'False')
plot_APD_vl(all_runs, 'F3', legend = 'False')

# Want to plot time (Age) versus APD, no simulation
plot_APD_time_ggplot <- function(data){
    data = as.data.table(data)
    data = data[, subject := .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]
    data = data[!is.na(apd)]

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))

    filename = paste0('figures/APD_versus_Time_by_Patient_for_all.png')

    # Create plot

    png(file=filename, width = 10, height = 5, units = 'in', res = 200)

    ggplot(data=data ,aes(y = apd, x=time, color=Sample)) + geom_point(size = 2.5) + facet_wrap(~Fragment) + geom_smooth(method="lm", se = FALSE, fullrange=TRUE) + ggtitle(paste0('APD versus Time by Patient--INFANT COHORT')) + labs(y="Average Pairwise Diversity (APD)", x= "Sampling Time (years)", color = "Patient")  + theme_bw() +  theme(text = element_text(size = 16)) + scale_fill_brewer() + ylim(0,0.04)

    dev.off()
}



# Want to plot time (Age) versus APD, no simulation
plot_APD_vl_ggplot <- function(data){
    data = as.data.table(data)
    data = data[, subject := .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]
    data = data[!is.na(apd)]
    data$vl = log10(data$VL)

    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$subject))

    filename = paste0('figures/APD_versus_VL_by_Patient_for_all.png')

    # Create plot

    png(file=filename, width = 10, height = 5, units = 'in', res = 200)

    ggplot(data=data ,aes(y = apd, x=vl)) + geom_point(aes(color=Sample), size = 2.5) + facet_wrap(~Fragment) + geom_smooth(method="lm", se = TRUE, fullrange=TRUE) + ggtitle(paste0('APD versus Viral Load--NEHER COHORT')) + labs(y="Average Pairwise Diversity (APD)", x= "log10 Viral Load", color = "Patient")  + stat_cor(method = "pearson", size = 4, label.x = 0, label.y = 0.0225) + stat_regline_equation(label.x = 0, label.y = 0.025, size = 4) + theme_bw() +  theme(text = element_text(size = 16)) + scale_fill_brewer()

    dev.off()
}

plot_APD_time_ggplot(all_runs)
plot_APD_vl_ggplot(all_runs)

# Want to plot time (Age) versus APD, no simulation
plot_APD_rate_vl_ggplot <- function(data){
    data = as.data.table(data)
    data = data[, subject := .GRP, by = Sample]
    data$time = data[["ActualTOI..year."]]
    data$apd = data[["APD_1"]]
    data_sp_vl = data[, sp_vl := log10(mean(VL)), by = Sample]
    data_rate = data.table()
    for (samp in unique(data$Sample)){
        for (frag in unique(data$Fragment)){
            data_temp = data[Sample == samp & Fragment == frag]
            if (nrow(data_temp)>0){
                data_temp$apd_rate = lm(apd ~ time, data = data_temp)[[1]][2]
                data_rate = rbind(data_rate, data_temp)
            }
        }
    }
    data_combined = unique(merge(data_rate, data_sp_vl)[,c(2,20, 37, 40,41)])
    

    par(mar=c(5,6,4,4)+.1)d
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data_combined$subject))

    filename = paste0('figures/APD_rate_versus_VL_by_Patient_for_all.png')

    # Create plot

    png(file=filename, width = 10, height = 5, units = 'in', res = 200)

    ggplot(data=data_combined ,aes(y = apd_rate, x=sp_vl)) + geom_point(aes(color=Sample), size = 2.5) + facet_wrap(~Fragment) + geom_smooth(method="lm", se = TRUE, fullrange=TRUE) + ggtitle(paste0('APD Rate versus Set Point Viral Load')) + labs(y="APD Rate of Change (diversity/year)", x= "log10 Set Point Viral Load", color = "Patient")  + stat_cor(method = "pearson", size = 4, label.x = 5.5, label.y = 0.0175) + stat_regline_equation(label.x = 5.5, label.y = 0.02, size = 4) + theme_bw() +  theme(text = element_text(size = 16)) + scale_fill_brewer()

    dev.off()
}

png(file=filename, width = 10, height = 5, units = 'in', res = 200)

    ggplot(data=data_combined ,aes(y = apd_rate, x=sp_vl)) + geom_point(aes(color=Sample), size = 2.5) + facet_wrap(~Fragment) + geom_smooth(method="lm", se = TRUE, fullrange=TRUE) + ggtitle(paste0('APD Rate versus Set Point Viral Load')) + labs(y="APD Rate of Change (diversity/year)", x= "log10 Set Point Viral Load", color = "Patient")  + stat_cor(method = "pearson", size = 4, label.x = 2, label.y = 0.005) + stat_regline_equation(label.x = 2, label.y = 0.0055, size = 4) + theme_bw() +  theme(text = element_text(size = 16)) + scale_fill_brewer()

    dev.off()
}