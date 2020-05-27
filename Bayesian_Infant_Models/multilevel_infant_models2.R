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


# This incorporates time error, but does not have subject specific time error (not accounting for different infection times)
multilevel_fragment_subject_var_slope_model_time_error <- map2stan(
    alist(
        # residuals
        infection_time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),

        # linear models
        population_mean <- total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],


        # Account for error in time measurements (convert to time since infection)
        time ~ dnorm(infection_time, time_error),

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 

        # fixed priors
        time_error ~ dunif(0, 0.75),
        subject_time_mean ~ dnorm(0, 0.1), 
        subject_time_sd ~ dcauchy(0,0.1),
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, start = list(infection_time = indexed_infant_data_cleaned_subset$time), WAIC = FALSE, iter = 20000, warmup = 5000, chains = 4, cores = 4, control = list(adapt_delta= 0.99)
)


# This incorporates time error for subject specific time error (accounting for different infection times)
multilevel_fragment_subject_var_slope_model_time_error2 <- map2stan(
    alist(
        # residuals
        infection_time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),

        # linear models
        population_mean <- total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],


        # Account for error in time measurements (convert to time since infection)
        time ~ dnorm(infection_time, time_error),

        time_error <- avg_time_error + subject_time_error[subject_index],

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        subject_time_error[subject_index] ~ dnorm(subject_time_mean, subject_time_sd),

        # fixed priors
        avg_time_error ~ dunif(0, 0.75),
        subject_time_mean ~ dnorm(0, 0.1), 
        subject_time_sd ~ dcauchy(0,0.1),
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, start = list(infection_time = indexed_infant_data_cleaned_subset$time +0.2), init_r = 0.5, WAIC = FALSE, iter = 30000, warmup = 10000, chains = 3, cores = 3, control = list(max_treedepth = 20, adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_fragment_subject_var_slope_model_time_error2", type = "time_error")


## Same time error method as above, but now allowing for varying intercepts

multilevel_fragment_subject_var_slope_model_time_error2_var_int <- map2stan(
    alist(
        # residuals
        infection_time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),

        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],


        # Account for error in time measurements (convert to time since infection)
        time ~ dnorm(infection_time, time_error),

        time_error <- avg_time_error + subject_time_error[subject_index],

        # adaptive priors
        subject_slope[subject_index] ~ dnorm(subject_slope_mean, subject_slope_sd), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        subject_time_error[subject_index] ~ dnorm(subject_time_mean, subject_time_sd),

        # fixed priors
        total_intercept ~ dnorm(0, 0.1),
        avg_time_error ~ dunif(0, 0.75),
        subject_time_mean ~ dnorm(0, 0.1), 
        subject_time_sd ~ dcauchy(0,0.1),
        subject_slope_mean ~ dnorm(0,1),
        subject_slope_sd ~ dcauchy(0,20),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, start = list(infection_time = indexed_infant_data_cleaned_subset$time +0.2), init_r = 0.5, WAIC = FALSE, iter = 30000, warmup = 10000, chains = 3, cores = 3, control = list(max_treedepth = 20, adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_fragment_subject_var_slope_model_time_error2_var_int", type = "time_error_int")



multilevel_fragment_subject_var_slope_model_time_error2_var_int_multivar_norm <- map2stan(
    alist(
        # residuals
        infection_time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,0.5),

        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_slope <- (population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index]) & T[0,],


        # Account for error in time measurements (convert to time since infection)
        time ~ dnorm(infection_time, time_error),

        time_error <- avg_time_error + subject_time_error[subject_index],

        # adaptive priors
        c(subject_slope, subject_time_error)[subject_index] ~ dmvnorm2(c(subject_slope_mean, subject_time_mean), subject_sd, subject_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 

        # fixed priors
        total_intercept ~ dnorm(0, 0.1),
        avg_time_error ~ dunif(0, 0.75),
        subject_time_mean ~ dnorm(0, 0.1), 
        subject_slope_mean ~ dnorm(0,1),
        subject_sd ~ dcauchy(0,20),
        subject_correlation_matrix ~ dlkjcorr(2),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,1),
        fragment_slope_sd ~ dcauchy(0,5)
    ), 
    data = indexed_infant_data_cleaned_subset, start = list(infection_time = indexed_infant_data_cleaned_subset$time +0.2), init_r = 0.5, WAIC = FALSE, iter = 30000, warmup = 10000, chains = 3, cores = 3, control = list(max_treedepth = 20, adapt_delta= 0.99)
)


map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,1),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_intercept <- population_avg_intercept + subject_intercept[subject_index],
        total_slope <- population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index],
        # adaptive priors
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnorm2(c(subject_intercept_change_mean, subject_slope_change_mean), subject_sd, subject_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        # fixed priors
        subject_intercept_change_mean ~ dnorm(0,0.2),
        subject_slope_change_mean ~ dnorm(0,0.5),
        subject_sd ~ dcauchy(0,20),
        subject_correlation_matrix ~ dlkjcorr(2),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,0.5),
        fragment_slope_sd ~ dcauchy(15,2)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 3, cores = 3, control = list(adapt_delta= 0.99)
)
