library(rethinking)
library(data.table)


data = read.csv("../_ignore/AllRunsAvg.csv")

preprocess <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA', .(mean(APD_1), mean(APD_5), mean(APD_10)), by = .(Sample, ActualTOI..year., Fragment, VL)]
    colnames(data_avg) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    return(data_avg)
}


# Trying to recapitualate OLS regressions...

# Fit model for lm style: 

fit_regression_model <- function(data, apd, fragment, type){
    data1 = preprocess(data)
    if (fragment == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (fragment == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (fragment == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    data_frag$apd = data_frag[[paste0('avg_apd', apd)]]
    if (type == 'LM'){
        model_fit_reg <- map(
        alist(
            time ~ dnorm(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(-0.75, 1), 
            b ~ dunif(0,150), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag
        ) 
    } else if (type == 'LAD') {
        model_fit_reg <- map(
        alist(
            time ~ dlaplace(mu, sigma), 
            mu <- a +b*apd, 
            a ~ dunif(-0.75, 1), 
            b ~ dunif(0,150), 
            sigma ~ dunif(0, 0.75)
            ), 
        data = data_frag
        ) 
    }
    return(model_fit_reg)
}


for (frag in c('F1', 'F2', 'F3')){
    for (ty in c('LM', "LAD")){
        assign(paste0(ty,'_model_bay1', frag), fit_regression_model(data, 1, frag, ty))
    }
}


plot_regression <- (data, type, apd_cut = 1, frag){
    model = paste0(type, '_model_bay', apd_cut, frag)
    data1 = preprocess(data)
    if (frag == "F1"){
        data_frag = data1[fragment == "F1"]
    } else if (frag == "F2"){
        data_frag = data1[fragment == "F2"]
    } else if (frag == "F3"){
        data_frag = data1[fragment == "F3"]
    }
    data_frag$apd = data_frag[[paste0('avg_apd', apd_cut)]]
    # Get a group of apd values in a sequence
    apd.seq <- seq(from = 0, to = 0.02, by = 0.001)
    # use the model to sample from the posterior for each apd value
    mu <- link(get(model), data = data.frame(apd = apd.seq))
    str(mu)

    plot(time ~ apd, data_frag, type = 'n')

    # Calculate the mean of the posterior distribution
    mu.mean <- apply(mu, 2, mean)

    # Calculate the standard deviation of the posterior distribution 
    #(contains 90% lower and upper bounds for each apd value)
    mu.HPDI <- apply(mu, 2, HPDI, prob = 0.9)

    # Plot the time versus apd for the data
    plot(time ~ apd, data = data_frag, col = col.alpha(rangi2, 0.5), main=paste0("APD", apd_cut, " versus Actual Time Since Infection"))

    # Add a line based on the mu's sampled for each apd value (i.e. the regression line)
    lines(apd.seq, mu.mean)

    # Add the confidence interval   
    shade(mu.HPDI, apd.seq)
    a = as.numeric(coef(get(model))['a'])
    b = as.numeric(coef(get(model))['b'])
    text(0.015, 0.25, paste0("y = ", as.numeric(round(coef(get(model))['a'], 4))," + ", round(as.numeric(coef(get(model))['b']), 4), "*x"))
}







