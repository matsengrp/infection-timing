library(data.table)
library(tidyverse)
library(rstan)
library(foreach)
library(cowplot)

TIME_CORRECTION_TYPE <<- 'beta' 
NCPU <<- 2
PREPROCESS_DATA <<- TRUE

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plotting_functions.R'))

infant_data = configure_data(TRAINING_INFANT_DATA_PATH)
data = as.data.table(infant_data)


# Measure relationship between APD and time over all patients and sequence regions

regression = lm(apd ~ observed_time_since_infection, data = data)
confint(regression, 'observed_time_since_infection', level=0.95)

# Quantify variation between subjects
subject_regression = lm(apd ~ observed_time_since_infection:as.factor(subject_id), data = data)
slopes_by_subject = as.data.frame(subject_regression$coefficients)
colnames(slopes_by_subject) = 'coef'
slopes_by_subject$term = rownames(slopes_by_subject)
slopes_by_subject = as.data.table(slopes_by_subject)
slopes = slopes_by_subject[!(term == '(Intercept)')]
range(slopes$coef)
anova(subject_regression)

# Quantify variation between sequence regions
frag_regression = lm(apd ~ observed_time_since_infection:as.factor(fragment), data = data)
slopes_by_frag = as.data.frame(frag_regression$coefficients)
colnames(slopes_by_frag) = 'coef'
slopes_by_frag$term = rownames(slopes_by_frag)
slopes_by_frag = as.data.table(slopes_by_frag)
frag_slopes = slopes_by_frag[!(term == '(Intercept)')]
range(frag_slopes$coef)
anova(frag_regression)

# Quantify variation between early/late
time_regression = lm(apd ~ observed_time_since_infection:as.factor(infection_status), data = data)
anova(time_regression)

sim_data = data[, min(observed_time_since_infection), by = subject_id] 
setnames(sim_data, 'V1', 'min_time')
sim_data = merge(sim_data, data[, max(observed_time_since_infection), by = subject_id])
setnames(sim_data, 'V1', 'max_time')
sim_data = sim_data[, seq(min_time, max_time, by = 0.01), by = subject_id]
setnames(sim_data, 'V1', 'observed_time_since_infection')
sim_data = sim_data[, predicted_apd := predict(subject_regression, sim_data)]

plot_apd_time_all(data, sim_data)
