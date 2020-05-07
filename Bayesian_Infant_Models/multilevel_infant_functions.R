library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

d = read.csv("../_ignore/AllRunsAvg.csv")

preprocess_noavg <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA']
    data_avg2 = data_avg %>% select(Sample, ActualTOI..year., Fragment, VL, APD_1, APD_5, APD_10)
    colnames(data_avg2) = c("sample", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    data_avg2$time = data_avg2$time + 0.125
    return(data_avg2)
}

dp = preprocess_noavg(d)
apd = 1
dp$apd = dp[[paste0('avg_apd', apd)]]
dp = dp[order(sample, fragment)]
sample_index = c(seq(1, length(unique(dp$sample))))
names(sample_index) = unique(dp$sample)
frag_index = c(seq(1, length(unique(dp$fragment))))
names(frag_index) = unique(dp$fragment)
indexed_dp = data.frame()
for (samp in unique(dp$sample)){
    dp1 = NULL
    dp1 = dp[with(dp, sample == samp)]
    dp1$sample_index = sample_index[[samp]]
    for (frag in unique(dp$fragment)){
        dp2 = NULL
        dp2 = dp1[with(dp1, fragment == frag)]
        dp2$frag_index = frag_index[[frag]]
        indexed_dp = rbind(indexed_dp, dp2)
    }   
}
#fit dummy variable model (for patient), but not multilevel

m12.1 <- map2stan(
    alist(
        time ~ dnorm(mu, sigma), 
        mu <- a_sample[sample_index] + b*apd,
        a_sample[sample_index] ~ dnorm(0, 0.25), 
        b ~ dunif(0, 150), 
        sigma ~ dcauchy(0,1)
    ), 
    data = indexed_dp
)

# inspect estimates (notice a different intercept for each tank)
precis(m12.1, depth = 2)

# plot distribution for visualization
x = seq(-0.75, 0.75, 0.01)
yfit<-dnorm(x,mean=-.375,sd=0.2) 

plot(x, yfit)
lines(x, yfit, col="blue", lwd=2)

# Now, make multilevel model with only patient level

m12.2 <- map2stan(
    alist(
        time ~ dnorm(mu, sigma), 
        mu <- a + a_sample[sample_index] + b*apd,
        a ~ dnorm(-0.125, 0.2),
        b ~ dunif(0, 150), 
        sigma ~ dcauchy(0,1),
        a_sample[sample_index] ~ dnorm(alph_sample, sigma_sample), 
        alph_sample ~ dnorm(0,0.2),
        sigma_sample ~ dcauchy(0,0.25)
    ), 
    data = indexed_dp, iter = 4000, chains = 4
)

# plot trace plot and posterior distribution to make sure the effective number of samples and Rhat values look alright
plot(m12.4)
precis(m12.4, depth = 2)

# Note that a_actor are DEVIATIONS from a, so here, compute the total intercept for each patient
post = extract.samples(m12.2)
total_a_sample = sapply(1:11, function(sample_index) post$a + post$a_sample[,sample_index])
round(apply(total_a_sample, 2, mean), 2)


# Add fragments as another level

m12.3 <- map2stan(
    alist(
        time ~ dnorm(mu, sigma), 
        mu <- a + a_sample[sample_index] + a_frag[frag_index] + b*apd,
        a ~ dnorm(-0.375, 0.2),
        b ~ dunif(0, 150), 
        sigma ~ dcauchy(0,1),
        a_sample[sample_index] ~ dnorm(alph_sample, sigma_sample), 
        alph_sample ~ dnorm(0,0.2),
        sigma_sample ~ dcauchy(0,0.25),
        a_frag[frag_index] ~ dnorm(alph_frag, sigma_frag), 
        alph_frag ~ dnorm(0,0.2),
        sigma_frag ~ dcauchy(0,0.25)
    ), 
    data = indexed_dp, iter = 4000, chains = 4
)

precis(m12.3, depth = 2)
plot(precis(m12.3, depth = 2))
plot(m12.3)


# compare sigma_actor and sigma_block
post <- extract.samples(m12.3)
dens(post$sigma_frag, xlab = 'sigma', xlim = c(0,4))
dens(post$sigma_sample, col = rangi2, lwd = 2, add = TRUE)
text(2, 0.85, 'sample', col = rangi2)
text(0.75, 2, 'frag')

# compare models
compare(m12.3, m12.2)

# predictions without link
post <- extract.samples(m12.3)
str(post)

dens(post$a_sample[,1])

## PREDICTIONS with training data
# create own link
p.link <- function(apd, frag_index, sample_index){
    time <- with(post, a + a_sample[,sample_index] + a_frag[,frag_index] + b * apd)
    return(time)
}

#compute predicitons
prediction_dp = data.frame()
for (samp in unique(indexed_dp$sample_index)){
    dp1 = NULL
    dp1 = indexed_dp[with(indexed_dp, sample_index == samp)]
    for (frag in unique(dp1$frag_index)){
        dp2 = NULL
        apd = NULL
        dp2 = dp1[with(dp1, frag_index == frag)]
        apd = dp2$apd
        pred.raw <- sapply(1:length(apd), function(i) p.link(apd[i], frag, samp))
        pred.p <- apply(pred.raw, 2, mean)
        pred.p.PI <- apply(pred.raw, 2, PI)
        dp2$pred_time = pred.p
        prediction_dp = rbind(prediction_dp, dp2)
    }   
}

plot_regression <- function(data){ 
    par(mar=c(5,6,4,1)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    plot(prediction_dp$apd, prediction_dp$time, col = c(factor(data$sample), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), ylab = 'Observed Time since infection', xlab = 'Average Pairwise Diversity', main = paste0('APD versus Observed Time Since Infection by Patient'), cex.main = 1, cex.lab = 1, cex.axis = 1)

    lines(prediction_dp$apd, prediction_dp$pred_time, lty = 2, lwd = 4)
    legend("topleft", legend=levels(factor(data$sample)), pch=16, col=unique(factor(data$sample)), ncol = 2, cex = 0.75)
}
plot_regression(prediction_dp)


plot_ETI_TI <- function(data){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    plot(prediction_dp$pred_time, prediction_dp$time, col = c(factor(data$sample), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), xlim = c(-.1, 2.5), ylab = 'Observed Time since infection', xlab = 'Estimated Time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient'), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    abline(lm(data$time ~ data$pred_time), col = 'red', lty = 2, lwd = 4)
    abline(v=0.125, col = c('blue', alpha = 0.5))
    abline(h = 0.125, col = c('blue', alpha = 0.5))
    legend("topleft", legend=levels(factor(data$sample)), pch=16, col=unique(factor(data$sample)), ncol = 2, cex = 0.75)
}
plot_ETI_TI(prediction_dp)

plot_ETI_TI_by_sample <- function(data){ 
    par(mar=c(5,6,4,4)+.1)
    palette(brewer.pal(n = 11, name = 'Set3'))
    col = setNames(palette(), levels(data$sample))
    plot(prediction_dp$pred_time, prediction_dp$time, col = c(factor(data$sample), alpha = 0.6), pch=19, cex = 1.25, ylim = c(0, 2), xlim = c(-1, 2.5), ylab = 'Observed Time since infection', xlab = 'Estimated Time since infection', main = paste0('Estimated versus Observed Time Since Infection by Patient'), cex.main = 1, cex.lab = 1, cex.axis = 1, panel.first = grid())
    for (samp in unique(data$sample)){
        data_s = data[with(data, sample == samp)]
        abline(lm(data_s$time ~ data_s$pred_time), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
    }
    abline(v=0.125, col = c('blue', alpha = 0.5))
    abline(h = 0.125, col = c('blue', alpha = 0.5))
    abline(a = 0, b = 1, lty = 2, lwd = 4)
    legend("topleft", legend=levels(factor(data$sample)), pch=16, col=unique(factor(data$sample)), ncol = 2, cex = 0.75)
}

plot_ETI_TI_by_sample(prediction_dp)

