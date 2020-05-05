library(rethinking)
library(data.table)
library(dplyr)
library(RColorBrewer)

for (frag in c('F1', 'F2', 'F3')){
    for (ty in c('LM', "LAD")){
        assign(paste0(ty,'_model_bay1', frag), readRDS(paste0(ty, "_model_bay1", frag, "_noAVG.rds")))
        assign(paste0(ty,'_model_bayHMC1', frag), readRDS(paste0(ty, "_model_bayHMC1", frag, ".rds")))
    }
}


preprocess <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10)), by = .(Sample, ActualTOI..year., Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    return(data_avg)
}

preprocess_noavg <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA']
    data_avg2 = data_avg %>% select(Sample, ActualTOI..year., Fragment, VL, APD_1, APD_5, APD_10)
    colnames(data_avg2) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    return(data_avg2)
}

plot_regression <- function(data, type, apd_cut, frag, method){ 
    data1 = preprocess_noavg(data)
    if (frag == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (frag == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (frag == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    data_frag$apd = data_frag[[paste0('avg_apd', apd_cut)]]
    # Get a group of apd values in a sequence
    apd.seq <- seq(from = 0, to = 0.025, by = 0.001)
    # use the model to sample from the posterior for each apd value
    if (method == 'quad'){
        mu <- link(get(paste0(type, '_model_bay', apd_cut, frag)), data = data.frame(apd = apd.seq))
        method_name = ' Normal Optimization'
    } else if (method == 'HMC'){
        mu <- link(get(paste0(type, '_model_bayHMC', apd_cut, frag)), data = data.frame(apd = apd.seq))
        method_name = ' Hamiltonian Monte Carlo Optimization'
    }
    # Calculate the mean of the posterior distribution
    mu.mean <- apply(mu, 2, mean)

    # Calculate the standard deviation of the posterior distribution 
    #(contains 90% lower and upper bounds for each apd value)
    mu.HPDI <- apply(mu, 2, HPDI, prob = 0.9)

    # Plot the time versus apd for the data
    par(mar=c(5,6,4,1)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    plot(time ~ apd, data = data_frag, col = factor(data_frag$sample), pch=19, cex = 1.25, main=paste0("APD", apd_cut, " versus Actual Time Since Infection \n",  "Fragment ", frag, " ", type, " Regression \n", method_name), cex.main = 1, cex.lab = 1, cex.axis = 1)
    legend("topleft", legend=levels(factor(data_frag$sample)), pch=16, col=unique(factor(data_frag$sample)), ncol = 2, cex = 0.75)
    # Add a line based on the mu's sampled for each apd value (i.e. the regression line)
    lines(apd.seq, mu.mean, lwd = 4)

    # Add the confidence interval   
    shade(mu.HPDI, apd.seq)
    if (method == 'quad'){
        a = as.numeric(coef(get(paste0(type, '_model_bay', apd_cut, frag)))['a'])
        b = as.numeric(coef(get(paste0(type, '_model_bay', apd_cut, frag)))['b'])
        text(0.015, 0.25, paste0("y = ", as.numeric(round(coef(get(paste0(type, '_model_bay', apd_cut, frag)))['a'], 4))," + ", round(as.numeric(coef(get(paste0(type, '_model_bay', apd_cut, frag)))['b']), 4), "*x"), cex = 1)
    } else if (method == 'HMC'){
        a = as.numeric(coef(get(paste0(type, '_model_bayHMC', apd_cut, frag)))['a'])
        b = as.numeric(coef(get(paste0(type, '_model_bayHMC', apd_cut, frag)))['b'])
            text(0.015, 0.25, paste0("y = ", as.numeric(round(coef(get(paste0(type, '_model_bayHMC', apd_cut, frag)))['a'], 4))," + ", round(as.numeric(coef(get(paste0(type, '_model_bayHMC', apd_cut, frag)))['b']), 4), "*x"), cex = 1)
    }
}

plot_ETI_TI <- function(data, type, apd_cut, frag, method){ 
    data1 = preprocess_noavg(data)
    if (frag == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (frag == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (frag == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    data_frag$apd = data_frag[[paste0('avg_apd', apd_cut)]]
    type_name = paste0(' Fragment', frag, ' ', type, ' Regression')
    # use the model to sample from the posterior for each apd value
    if (method == 'quad'){
        mu <- link(get(paste0(type, '_model_bay', apd_cut, frag)))
        time.sim <- sim(get(paste0(type, '_model_bay', apd_cut, frag)), n= 1e4)
        method_name = ' Normal Optimization'
    } else if (method == 'HMC'){
        mu <- link(get(paste0(type, '_model_bayHMC', apd_cut, frag)))
        time.sim <- sim(get(paste0(type, '_model_bayHMC', apd_cut, frag)), n= 1e4)
        method_name = ' Hamiltonian Monte Carlo Optimization'
    }
    # summarize samples across cases
    mu.mean <- apply(mu, 2, mean)
    mu.PI <- apply(mu, 2, PI)

    # simulate observations
    time.PI <- apply(time.sim, 2, PI)
    par(mar=c(5,6,4,1)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    plot(mu.mean ~ data_frag$time, col = c(factor(data_frag$sample), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), xlab = 'Observed Time since infection', ylab = 'Estimated Time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient \n', type_name, "\n", method_name), cex.main = 1, cex.lab = 1, cex.axis = 1)
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    for (i in 1:nrow(data_frag)){
        lines(rep(data_frag$time[i], 2), c(mu.PI[1,i], mu.PI[2, i]), col = alpha("black", 0.4), lwd = 2)
    }
    legend("topleft", legend=levels(factor(data_frag$sample)), pch=16, col=unique(factor(data_frag$sample)), ncol = 2, cex = 0.75)
}

model_eval <- function(model1, model2){
    models = compare(get(model1), get(model2))
    plot(models, SE = TRUE, dSE = TRUE)
    plot(coeftab(get(model1), get(model2)))
}
    

