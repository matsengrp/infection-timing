library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("multilevel_infant_functions.R")

# No divergent iterations here!
multilevel_fragment_subject_var_slope_model_new2 <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],

        # Intercept 
        total_intercept ~ dnorm(0, 0.5) & T[-0.75, 0],

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 

        # fixed priors
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 4, cores = 4, control = list(adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_fragment_subject_var_slope_model_new", type = "fragment_subject_varying_slopes")



multilevel_subject_var_int_fragment_subject_var_slope_model_new2 <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]), 
        total_intercept <- (population_avg_int + subject_intercept[subject_index]), 

        # Intercept 

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd),
        subject_intercept[subject_index] ~ dnorm(subject_intercept_mean, subject_intercept_sd),
        # fixed priors
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        subject_intercept_mean ~ dnorm(0,1),
        subject_intercept_sd ~ dcauchy(0,0.25),
        population_avg_slope ~ dunif(0, 150),
        population_avg_int  ~ dnorm(0,0.5),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 4, cores = 4, control = list(adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_fragment_subject_var_slope_model_new", type = "fragment_subject_varying_slopes")


# No divergent iterations here!
multilevel_fragment_subject_var_slope_model_new_no_int <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),
        # linear models
        population_mean <- total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 

        # fixed priors
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 4, cores = 4, control = list(adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_fragment_subject_var_slope_model_new_no_int", type = "fragment_subject_varying_slopes_no_int")
