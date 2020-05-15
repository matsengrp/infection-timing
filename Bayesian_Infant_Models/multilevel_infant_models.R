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


# Now, make multilevel model with only patient level ( this model does not have divergence problems! Woo!)

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


# Add fragments as another level for varying intercepts ( this model gives divergent solutions, but I am not going to pursue it, because it does not make sense as a model to have fragment dictate intercept, fragment should, instead, influence slope...)

#multilevel_subject_fragment_var_int_model <- map2stan(
#    alist(
#        time ~ dnorm(population_mean, population_sd), 
#        population_mean <- population_avg_intercept + subject_intercept[subject_index] + #fragment_intercept[frag_index] + population_avg_slope*apd,
#        population_avg_intercept ~ dunif(-0.75, 0),
#        population_avg_slope ~ dunif(0, 150), 
#        population_sd ~ dcauchy(0,1),
#        subject_intercept[subject_index] ~ dnorm(subject_mean, subject_sd), 
#        subject_mean ~ dnorm(0,0.2),
#        subject_sd ~ dcauchy(0,0.25),
#        fragment_intercept[frag_index] ~ dnorm(fragment_mean, fragment_sd), 
#        fragment_mean ~ dnorm(-50,50),
#        fragment_sd ~ dcauchy(0,1)
#    ), 
#    data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 4, #cores = 4, refresh = 20000, control = list(adapt_delta= 0.99)
#)
#
#save_models_plots_stats("multilevel_subject_fragment_var_int_model", type = #"subject_fragment_varying_intercepts")


# Add varying slope for fragments and keep varying intercepts for subjects as another level (no divergences here! Woo!)

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
        fragment_mean ~ dnorm(0,0.5),
        fragment_sd ~ dcauchy(15,2)
    ), 
   data = indexed_infant_data_cleaned_subset, iter = 20000, warmup = 5000, chains = 4, cores = 4, refresh = 20000, control = list(adapt_delta= 0.99)
)

save_models_plots_stats("multilevel_subject_var_int_fragment_var_slope_model", type = "subject_varying_intercepts_fragment_varying_slopes")


# Add varying slope for subjects and keep varying intercepts for subjects only as another level (this one has many divergences (500 or so))

#multilevel_subject_var_int_subject_var_slope_model <- #map2stan(
#    alist(
#        # residuals
#        time ~ dnorm(population_mean, population_sd), 
#        population_sd ~ dcauchy(0,1),
#        # linear models
#        population_mean <- total_intercept + total_slope*apd,
#        total_intercept <- population_avg_intercept + #subject_intercept[subject_index],
#        total_slope <- population_avg_slope + subject_slope#[subject_index],
#        # adaptive priors
#        c(subject_intercept, subject_slope)[subject_index] ~ #dmvnorm2(c(subject_intercept_change_mean, #subject_slope_change_mean), subject_sd, #subject_correlation_matrix), 
#        # fixed priors
#        # for adaptive prior
#        subject_intercept_change_mean ~ dnorm(0,0.2),
#        subject_slope_change_mean ~ dnorm(0,0.5),
#        subject_sd ~ dcauchy(0,20),
#        subject_correlation_matrix ~ dlkjcorr(1),
#        population_avg_intercept ~ dunif(-0.75, 0),
#        population_avg_slope ~ dunif(0, 150)
#    ), 
#   data = indexed_infant_data_cleaned_subset, iter = 20000, #warmup = 5000, chains = 3, cores = 3, control = list#(adapt_delta= 0.99)
#)
#indexed_infant_data_cleaned_subset[!(subject_index %in% c(6, #7, 9)),]
#save_models_plots_stats#("multilevel_subject_var_int_subject_var_slope_model", type #= "subject_varying_intercepts_slopes")

# Add varying slope for fragments and subjects and keep varying intercepts for subjects only as another level

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

save_models_plots_stats("multilevel_subject_var_int_fragment_subject_var_slope_model", type = "subject_varying_intercepts_fragment_subject_varying_slopes")

# Try non centered parameterization
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
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnormNC(subject_sd, subject_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        # fixed priors
        subject_sd ~ dcauchy(0,20),
        subject_correlation_matrix ~ dlkjcorr(2),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,0.5),
        fragment_slope_sd ~ dcauchy(15,2)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 30000, warmup = 5000, chains = 5, cores = 5, control = list(adapt_delta= 0.99)
)



# Now add another level for the multiple runs of apd measures to varying slope. So, in total, we will have subject index varying intercept and fragment index, subject index, and multiple run index all varying slope. 

multilevel_subject_var_int_fragment_subject_run_var_slope_model <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,1),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_intercept <- population_avg_intercept + subject_intercept[subject_index],
        total_slope <- population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index] + run_slope[run],
        # adaptive priors
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnorm2(c(subject_intercept_change_mean, subject_slope_change_mean), subject_sd, subject_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        run_slope[run] ~ dnorm(run_slope_mean, run_slope_sd), 
        # fixed priors
        # for adaptive prior
        subject_intercept_change_mean ~ dunif(-0.25, 0.25),
        subject_slope_change_mean ~ dunif(-50, 50),
        subject_sd ~ dcauchy(0, 0.25),
        subject_correlation_matrix ~ dlkjcorr(1),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dunif(-50,50),
        fragment_slope_sd ~ dcauchy(0,1),
        run_slope_mean ~ dunif(-50,50),
        run_slope_sd ~ dcauchy(0,1)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 10000, warmup = 2000, chains = 4, cores = 4, control = list(adapt_delta= 0.95)
)

save_models_plots_stats("multilevel_subject_var_int_fragment_subject_run_var_slope_model", type = "subject_varying_intercepts_fragment_subject_run_varying_slopes")

# Try non centered parameterization
multilevel_subject_var_int_fragment_subject_run_var_slope_model <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,1),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_intercept <- population_avg_intercept + subject_intercept[subject_index],
        total_slope <- population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index] + run_slope[run],
        # adaptive priors
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnormNC(subject_sd, subject_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        run_slope[run] ~ dnorm(run_slope_mean, run_slope_sd), 
        # fixed priors
        subject_sd ~ dcauchy(0,20),
        subject_correlation_matrix ~ dlkjcorr(2),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,0.5),
        fragment_slope_sd ~ dcauchy(15,2), 
        run_slope_mean ~ dnorm(0,0.5),
        run_slope_sd ~ dcauchy(18,2)
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 30000, warmup = 5000, chains = 5, cores = 5, control = list(adapt_delta= 0.99)
)

# Try non centered parameterization
multilevel_subject_run_var_int_fragment_subject_run_var_slope_model <- map2stan(
    alist(
        # residuals
        time ~ dnorm(population_mean, population_sd), 
        population_sd ~ dcauchy(0,1),
        # linear models
        population_mean <- total_intercept + total_slope*apd,
        total_intercept <- population_avg_intercept + subject_intercept[subject_index] + run_intercept[run],
        total_slope <- population_avg_slope + subject_slope[subject_index] + fragment_slope[frag_index] + run_slope[run],
        # adaptive priors
        c(subject_intercept, subject_slope)[subject_index] ~ dmvnormNC(subject_sd, subject_correlation_matrix),
        c(run_intercept, run_slope)[run_index] ~ dmvnormNC(run_sd, run_correlation_matrix), 
        fragment_slope[frag_index] ~ dnorm(fragment_slope_mean, fragment_slope_sd), 
        run_slope[run] ~ dnorm(run_slope_mean, run_slope_sd), 
        # fixed priors
        subject_sd ~ dcauchy(0,20),
        subject_correlation_matrix ~ dlkjcorr(2),
        run_sd ~ dcauchy(0,20),
        run_correlation_matrix ~ dlkjcorr(2),
        population_avg_intercept ~ dunif(-0.75, 0),
        population_avg_slope ~ dunif(0, 150),
        fragment_slope_mean ~ dnorm(0,0.5),
        fragment_slope_sd ~ dcauchy(15,2), 
    ), 
    data = indexed_infant_data_cleaned_subset, iter = 30000, warmup = 5000, chains = 5, cores = 5, control = list(adapt_delta= 0.99)
)

