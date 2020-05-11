library(rethinking)
library(rstan)
library(data.table)
library(dplyr)
library(RColorBrewer)

source("multilevel_infant_functions.R")

#fit dummy variable model (for patient), but not multilevel (THIS ISN'T WORKING FOR SOME REASON)

dummy_variable_subject_model <- map2stan(
    alist(
        time ~ dnorm(population_mean, population_sd), 
        population_mean <- subject_intercept[subject_index] + population_slope*apd,
        subject_intercept[subject_index] ~ dunif(-0.75, 0), 
        population_slope ~ dunif(0, 150), 
        population_sd ~ dcauchy(0,1)
    ), 
    data = indexed_infant_data_cleaned_subset2
)

# inspect estimates (notice a different intercept for each tank)
precis(dummy_variable_subject_model, depth = 2)


# Now, make multilevel model with only patient level

multilevel_subject_var_int_model <- map2stan(
    alist(
        time ~ dnorm(population_mean, population_sd), 
        population_mean <- population_avg_intercept + subject_intercept[subject_index] + population_avg_slope*apd,
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150), 
        population_sd ~ dcauchy(0,1),
        subject_intercept[subject_index] ~ dnorm(subject_mean, subject_sd),  
        subject_mean ~ dnorm(0,0.2),
        subject_sd ~ dcauchy(0,0.25)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 4000, chains = 4
)

save_models_plots_stats("multilevel_subject_var_int_model", type = "subject_varying_intercepts")


# Add fragments as another level for varying intercepts

multilevel_subject_fragment_var_int_model <- map2stan(
    alist(
        time ~ dnorm(population_mean, population_sd), 
        population_mean <- population_avg_intercept + subject_intercept[subject_index] + fragment_intercept[frag_index] + population_avg_slope*apd,
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150), 
        population_sd ~ dcauchy(0,1),
        subject_intercept[subject_index] ~ dnorm(subject_mean, subject_sd), 
        subject_mean ~ dnorm(0,0.2),
        subject_sd ~ dcauchy(0,0.25),
        fragment_intercept[frag_index] ~ dnorm(fragment_mean, fragment_sd), 
        fragment_mean ~ dnorm(0,0.2),
        fragment_sd ~ dcauchy(0,0.25)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 15000, warmup = 2000, chains = 4, cores = 4, control = list(adapt_delta= 0.95)
)

save_models_plots_stats("multilevel_subject_fragment_var_int_model", type = "subject_fragment_varying_intercepts")


# Add varying slope for fragments and keep varying intercepts for subjects as another level

multilevel_subject_var_int_fragment_var_slope_model <- map2stan(
    alist(
        time ~ dnorm(population_mean, population_sd), 
        population_mean <- population_avg_intercept + subject_intercept[subject_index] + (fragment_slope[frag_index] + population_avg_slope)*apd,
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150), 
        population_sd ~ dcauchy(0,1),
        subject_intercept[subject_index] ~ dnorm(subject_mean, subject_sd), 
        subject_mean ~ dnorm(0,0.2),
        subject_sd ~ dcauchy(0,0.25),
        fragment_slope[frag_index] ~ dnorm(fragment_mean, fragment_sd), 
        fragment_mean ~ dnorm(-50,50),
        fragment_sd ~ dcauchy(0,1)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 15000, warmup = 2000, chains = 4, cores = 4, control = list(adapt_delta= 0.95)
)

save_models_plots_stats("multilevel_subject_var_int_fragment_var_slope_model", type = "subject_varying_intercepts_fragment_varying_slopes")


# Add varying slope for subjects and keep varying intercepts for subjects only as another level

multilevel_subject_var_int_subject_var_slope_model <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,1),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_intercept <- population_avg_intercept + subject_intercept[subject_index],
        total_slope <- population_avg_slope + subject_slope[subject_index],
        # adaptive priors
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnorm2(c(subject_intercept_change_mean, subject_slope_change_mean), subject_sd, subject_correlation_matrix), 
        # fixed priors
        # for adaptive prior
        subject_intercept_change_mean ~ dunif(-0.25, 0.25),
        subject_slope_change_mean ~ dunif(-50, 50),
        subject_sd ~ dcauchy(0, 0.25),
        subject_correlation_matrix ~ dlkjcorr(1),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 10000, warmup = 2000, chains = 4, cores = 4, control = list(adapt_delta= 0.95)
)

save_models_plots_stats("multilevel_subject_var_int_subject_var_slope_model", type = "subject_varying_intercepts_slopes")

# Add varying slope for fragments and subjects and keep varying intercepts for subjects only as another level (THIS DOESNT WORK!!!)

multilevel_subject_var_int_fragment_subject_var_slope_model <- map2stan(
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
        # for adaptive prior
        subject_intercept_change_mean ~ dunif(-0.25, 0.25),
        subject_slope_change_mean ~ dunif(-50, 50),
        subject_sd ~ dcauchy(0, 0.25),
        subject_correlation_matrix ~ dlkjcorr(1),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dunif(-50,50),
        fragment_slope_sd ~ dcauchy(0,1)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 10000, warmup = 2000, chains = 4, cores = 4, control = list(adapt_delta= 0.95)
)

save_models_plots_stats("multilevel_subject_var_int_fragment_subject_var_slope_model", type = "subject_varying_intercepts_fragment_subject_varying_slopes")


