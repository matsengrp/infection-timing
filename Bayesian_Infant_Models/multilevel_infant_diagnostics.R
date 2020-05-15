library(rstan)
library(rethinking)
library(data.table)
library(dplyr)
library(RColorBrewer)

# need to load models
fit = multilevel_subject_fragment_var_int_model


if ( class(fit) %in% c("map2stan","ulam") ) fit <- fit@stanfit
    params_fit <- as.data.frame(extract(fit, permuted=FALSE))
    names(params_fit) <- gsub("chain:1.", "", names(params_fit), fixed = TRUE)
    names(params_fit) <- gsub("[", ".", names(params_fit), fixed = TRUE)
    names(params_fit) <- gsub("]", "", names(params_fit), fixed = TRUE)
    divergent <- rstan::get_sampler_params(fit, inc_warmup=FALSE)[[1]][,'divergent__']
    params_fit$iter <- 1:nrow(params_fit)

params_fit$divergent = divergent


div_params_fit <- params_fit[params_fit$divergent == 1,]
nondiv_params_fit <- params_fit[params_fit$divergent == 0,]

par(mar = c(4, 4, 0.5, 0.5))
plot(nondiv_params_fit$fragment_intercept.1, log(nondiv_params_fit$fragment_sd),
     xlab="fragment_intercept.1", ylab="fragment_sd",
     col="blue", pch=16, cex=0.8)
points(div_params_fit$fragment_intercept.1, log(div_params_fit$fragment_sd),
       col="green", pch=16, cex=0.8)

running_means_fit <- sapply(1:1800, function(n) mean(log(params_fit$fragment_sd)[1:(10*n)]))

plot(10*(1:1800), running_means_fit, col="blue", pch=16, cex=0.8,
    xlab="Iteration", ylab="MCMC mean of fragment_sd")
points(10*(1:1800), running_means_fit, col="blue", pch=16, cex=0.8)
abline(h=0, col="grey", lty="dashed", lwd=3)

par(mar = c(4, 4, 0.5, 0.5))
plot(params_fit$iter, log(params_fit$fragment_sd), col="blue", pch=16, cex=0.8,
     xlab="Iteration", ylab="fragment_sd")


### SIMULATE DATA
#
##simulate some data
#set.seed(20161110)
#N<-120 #sample size
#n_subject<-10 #number of subjects
#subject_id<-rep(1:n_subject,each=9) #index of subject
#n_fragments <- 3
#fragments = c(1:n_fragments)
#fragment_id<- rep(fragments, 9*4)
#apd = c(0.001, 0.004, 0.007)
#apds <- rep(apd, each = 3, 12)
#K<-3 #number of regression coefficients
##population-level regression coefficient
#gamma<-c(2,-1,3)
##standard deviation of the group-level coefficient
#tau<-c(0.3,2,1)
##standard deviation of individual observations
#sigma<-1
##group-level regression coefficients
#beta<-mapply(function(g,t) rnorm(J,g,t),g=gamma,t=tau) 
##the model matrix
#X<-model.matrix(~x+y,data=data.frame(x=runif(N,-2,2),y=runif(N,-2,2)))
#y<-vector(length = N)
#for(n in 1:N){
#  #simulate response data
#  y[n]<-rnorm(1,X[n,]%*%beta[id[n],],sigma)
#}



data.predicted = list(apd = c(0.001, 0.004, 0.007), subject_index = rep(1,3), frag_index = rep(1,3))
post <- extract.samples(multilevel_subject_var_int_fragment_var_slope_model)
post <- extract.samples(multilevel_subject_var_int_model)
subject_intercept_simulation = rnorm(7000, 0, post$subject_sd)
subject_intercept_simulation = matrix(subject_intercept_simulation, 1000, 7)
fragment_slope_simulation = rnorm(7000, 0, post$fragment_sd)
fragment_slope_simulation = matrix(fragment_slope_simulation, 1000, 7)

link.fit.frag = link(multilevel_subject_fragment_var_int_model, n = 1000, data = data.predicted, replace = list(subject_intercept = subject_intercept_simulation, fragment_slope = fragment_slope_simulation))

link.fit.frag = link(multilevel_subject_var_int_fragment_var_slope_model, n = 1000, data = data.predicted, replace = list(subject_intercept = subject_intercept_simulation, fragment_slope = fragment_slope_simulation))

link.fit = link(multilevel_subject_var_int_model, n = 1000, data = data.predicted, replace = list(subject_intercept = subject_intercept_simulation))

simulate.subject.frag <- function(i, apd){
    simulate_subject = rnorm(1,post$subject_mean[i], post$subject_sd[i])
    simulate_fragment = rnorm(1,post$fragment_mean[i], post$fragment_sd[i])
    time = post$population_avg_intercept[i] + simulate_subject + 
    (simulate_fragment + post$population_avg_slope[i]) * apd
    return(time)
}

simulate.subject <- function(i, apd){
    simulate_subject = rnorm(1,post$subject_mean[i], post$subject_sd[i])
    time = post$population_avg_intercept[i] + simulate_subject + post$population_avg_slope[i] * apd
    return(time)
}

# the following is great!!!!
apd = seq(0.0001, 0.03, by = 0.0003)
simulation = c()
for (j in apd){
    for (i in 1:50) simulation = c(simulation, simulate.subject.frag(i, j))
}
apd_sim = rep(apd, each = 50)
palette(brewer.pal(n = 11, name = 'Set3'))
plot(apd_sim,simulation, col = alpha("grey", 0.1), pch = 19, xlab = 'apd', ylab = 'time', xlim = c(0,0.03), ylim = c(0, max(simulation)+ 0.5))
points(indexed_infant_data_cleaned$apd, indexed_infant_data_cleaned$time, col = c(factor(indexed_infant_data_cleaned$subject), alpha = 0.6), pch=19, cex = 1.25)
for (samp in unique(indexed_infant_data_cleaned$subject)){
        data_s = indexed_infant_data_cleaned[with(indexed_infant_data_cleaned, subject == samp)]
        abline(lm(data_s$time ~ data_s$apd), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
    }
abline(lm(simulation~apd_sim), col="red", lw = 4) 

palette(brewer.pal(n = 11, name = 'Set3'))
col = setNames(palette(), levels(indexed_infant_data_cleaned$subject))
plot(indexed_infant_data_cleaned$apd, indexed_infant_data_cleaned$time, col = c(col, alpha = 0.6), pch = 19, xlab = 'apd', ylab = 'time', cex = 1.25)
for (samp in unique(indexed_infant_data_cleaned$subject)){
        data_s = indexed_infant_data_cleaned[with(indexed_infant_data_cleaned, subject == samp)]
        abline(lm(data_s$time ~ data_s$apd), col = c(col[[samp]], alpha = 0.5), lwd = 2.5)
    }
legend("topleft", legend=levels(factor(indexed_infant_data_cleaned$subject)), pch=16, col=unique(factor(indexed_infant_data_cleaned$subject)), ncol = 2, cex = 0.75)