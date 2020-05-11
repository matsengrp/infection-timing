library(rethinking)
library(data.table)
library(dplyr)


data = read.csv("../_ignore/AllRunsAvg.csv")

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



# Trying to recapitualate OLS regressions...

# Fit model for lm style: 

fit_regression_model <- function(data, apd, fragment, type, standardize){
    data1 = preprocess_noavg(data)
    if (fragment == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (fragment == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (fragment == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    if (standardize == 'FALSE'){
        data_frag$apd = data_frag[[paste0('avg_apd', apd)]]
    } else if (standardize == 'TRUE'){
        data_frag$apd = (data_frag[[paste0('avg_apd', apd)]]- mean(data_frag[[paste0('avg_apd', apd)]]))/sd(data_frag[[paste0('avg_apd', apd)]])
    }
    if (type == 'LM'){
        model_fit_reg <- map(
        alist(
            time ~ dnorm(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(0,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag, 
        method = 'Nelder-Mead', 
        control = list(maxit = 1e4)
        ) 
    } else if (type == 'LAD') {
        model_fit_reg <- map(
        alist(
            time ~ dlaplace(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(0,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag, 
        method = 'Nelder-Mead', 
        control = list(maxit = 1e4)
        ) 
    }
    return(model_fit_reg)
}

fit_regression_model_HMC <- function(data, apd, fragment, type, standardize){
    require(rstan)
    data1 = preprocess_noavg(data)
    if (fragment == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (fragment == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (fragment == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    data_frag$apd = data_frag[[paste0('avg_apd', apd)]]
    data_frag = data_frag %>% select(time, apd)
    if (type == 'LM'){
        model_fit_reg <- map2stan(
        alist(
            time ~ dnorm(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(30,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag
        ) 
    } else if (type == 'LAD') {
        model_fit_reg <- map(
        alist(
            time ~ dlaplace(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(30,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag
        ) 
    }
    return(model_fit_reg)
}

for (frag in c('F1', 'F2', 'F3')){
    for (ty in c('LM', "LAD")){
        assign(paste0(ty,'_model_bay1', frag), fit_regression_model(data, 1, frag, ty, standardize = FALSE))
        print(precis(get(paste0(ty, "_model_bay1", frag)), corr = TRUE))
        print(cov2cor(vcov(get(paste0(ty, "_model_bay1", frag)))))
        print("done")
        saveRDS(get(paste0(ty, "_model_bay1", frag)), paste0(ty, "_model_bay1", frag, "_noAVG.rds"))
    }
}

for (frag in c('F1', 'F2', 'F3')){
    for (ty in c('LM', "LAD")){
        assign(paste0(ty,'_model_bayHMC1', frag), fit_regression_model_HMC(data, 1, frag, ty, standardize = FALSE))
        print(precis(get(paste0(ty, "_model_bayHMC1", frag)), corr = TRUE))
        print(cov2cor(vcov(get(paste0(ty, "_model_bayHMC1", frag)))))
        print("done")
        saveRDS(get(paste0(ty, "_model_bayHMC1", frag)), paste0(ty, "_model_bayHMC1", frag, ".rds"))
    }
}

fit_regression_model_frag_apd <- function(data, apd, fragment, type, standardize){
    data1 = preprocess_noavg(data)
    if (fragment == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (fragment == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (fragment == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    if (standardize == 'FALSE'){
        data_frag$apd = data_frag[[paste0('avg_apd', apd)]]
    } else if (standardize == 'TRUE'){
        data_frag$apd = (data_frag[[paste0('avg_apd', apd)]]- mean(data_frag[[paste0('avg_apd', apd)]]))/sd(data_frag[[paste0('avg_apd', apd)]])
    }
    if (type == 'LM'){
        model_fit_reg <- map(
        alist(
            time ~ dnorm(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(0,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag, 
        method = 'Nelder-Mead', 
        control = list(maxit = 1e4)
        ) 
    } else if (type == 'LAD') {
        model_fit_reg <- map(
        alist(
            time ~ dlaplace(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(0, 0.75), 
            b ~ dunif(0,100), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag, 
        method = 'Nelder-Mead', 
        control = list(maxit = 1e4)
        ) 
    }
    return(model_fit_reg)
}






