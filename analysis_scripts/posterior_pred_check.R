library(rstan)
library(rstanarm)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)

TIME_CORRECTION_TYPE <<- 'beta'
NCPU <<- 2
PREPROCESS_DATA <<-  TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_cross_validation_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/calculation_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = infant_data
data$number = seq(1, length(data$is_post))
data = as.data.table(data)
data$fragment_id = data$fragment

model = load_model_fit()
post = rstan::extract(model)
transformed_posteriors = transform_posterior_matrix_to_dataframe(data, post$observed_time_since_infection_rep, observed_time = TRUE)
raw = unique(transformed_posteriors[, -c('iteration', 'time_since_infection_posterior_draw')])

for (stat in c('median', 'mean', 'sd', 'max', 'min')){
    print(calculate_posterior_pred_check(get(stat), transformed_posteriors, raw))
    plot_posterior_pred_check(get(stat), stat, transformed_posteriors, raw)
}

# adapt check for 'greater than'
iterations = length(unique(transformed_posteriors$iteration))
obs_length = length(raw$observed_time_since_infection)
post_stat = transformed_posteriors[time_since_infection_posterior_draw >= observed_time_since_infection, .N/obs_length, by = iteration]
stat = 0.5
p = nrow(post_stat[V1 >= stat])/iterations
greater_than_check = list(sim_stat = post_stat, stat = stat, p = p)
print(greater_than_check)

plot_posterior_pred_check('prop_greater-than', 'prop_greater-than', transformed_posteriors, raw, check = greater_than_check)

######## check individual observation coverage ##################

transformed_posteriors[, q0.055 := quantile(time_since_infection_posterior_draw, 0.055), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]

transformed_posteriors[, q0.945 := quantile(time_since_infection_posterior_draw, 0.945), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]

transformed_posteriors[, q0.255 := quantile(time_since_infection_posterior_draw, 0.255), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]

transformed_posteriors[, q0.745 := quantile(time_since_infection_posterior_draw, 0.745), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]

transformed_posteriors[, q0.405 := quantile(time_since_infection_posterior_draw, 0.405), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]

transformed_posteriors[, q0.595 := quantile(time_since_infection_posterior_draw, 0.595), by = .(apd, subject_id, fragment_id, observed_time_since_infection)]


cols = c(colnames(raw), 'q0.055', 'q0.945', 'q0.255', 'q0.745', 'q0.405', 'q0.595')

intervals = unique(transformed_posteriors[, ..cols])

# how many observations are covered by posterior simulation (within 0.89 interval)
intervals[observed_time_since_infection > q0.055 & observed_time_since_infection < q0.945, covered_89 := TRUE]
intervals[observed_time_since_infection <= q0.055 | observed_time_since_infection >= q0.945, covered_89 := FALSE]

print(intervals[, .N, by = covered_89])

# how many observations are covered by posterior simulation (within 0.49 interval)
intervals[observed_time_since_infection > q0.255 & observed_time_since_infection < q0.745, covered_49 := TRUE]
intervals[observed_time_since_infection <= q0.255 | observed_time_since_infection >= q0.745, covered_49 := FALSE]

print(intervals[, .N, by = covered_49])

# how many observations are covered by posterior simulation (within 0.19 interval)
intervals[observed_time_since_infection > q0.405 & observed_time_since_infection < q0.595, covered_19 := TRUE]
intervals[observed_time_since_infection <= q0.405 | observed_time_since_infection >= q0.595, covered_19 := FALSE]

print(intervals[, .N, by = covered_19])
