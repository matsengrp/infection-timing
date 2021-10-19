library(data.table)

source('../scripts/model_cross_validation_functions.R')

path = '../_ignore/neher_webtool_results.tsv'
adult_model_estimates = fread(path) 

# get time difference
adult_model_estimates[, time_abs_difference := abs(adult_model_time_since_infection_estimates - year_visit)]

# get mean absolute difference
print(mean(adult_model_estimates$time_abs_difference))
