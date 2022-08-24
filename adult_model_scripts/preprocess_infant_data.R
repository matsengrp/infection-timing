library(data.table)
library(tidyverse)
library(rstan)
library(foreach)

args = commandArgs(trailingOnly = TRUE)

DATA_PATH <<- args[1] 
TIME_KNOWN <<- args[2]
NCPU <<- 2
TIME_CORRECTION_TYPE <<- 'beta'
if (isTRUE(TIME_KNOWN)){
    PREPROCESS_DATA <<- TRUE
} else {
    PREPROCESS_DATA <<- FALSE 
}

source('config/config.R')
source(paste0(PROJECT_PATH, '/config/file_paths.R'))
source(paste0(PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(PROJECT_PATH, '/scripts/model_prediction_evaluation_functions.R'))

data = check_data(DATA_PATH, TIME_KNOWN)
# remove NA cases
data = data[!(is.na(pass2_APD))]
if (isTRUE(PREPROCESS_DATA)){
    data = preprocess_data(data)
}

if (isTRUE(TIME_KNOWN)){
    temp_cols = c('ptnum', 'month_visit', 'pass2_APD', 'Fragment', 'replicate', 'incat_hiv', 'inftimemonths', 'HXB2nt_start', 'HXB2nt_end')
} else {
    temp_cols = c('ptnum','pass2_APD', 'Fragment', 'replicate', 'HXB2nt_start', 'HXB2nt_end')
}

important_cols = temp_cols[!(temp_cols %in% c('replicate', 'pass2_APD', 'HXB2nt_start', 'HXB2nt_end'))]

subset = data[, ..temp_cols]

# average APD across replicates
average_apd = subset[, mean(pass2_APD), by = important_cols]
setnames(average_apd, 'V1', 'average_APD')
average_start = subset[, mean(HXB2nt_start), by = important_cols]
setnames(average_start, 'V1', 'HXB2nt_start')
average_end = subset[, mean(HXB2nt_end), by = important_cols]
setnames(average_end, 'V1', 'HXB2nt_end')

average_subset = merge(average_apd, average_start, by = important_cols)
average_subset = merge(average_subset, average_end, by = important_cols)

if (isTRUE(TIME_KNOWN)){
    # fill in missing infection times
    known_times = unique(subset[, c('ptnum', 'incat_hiv', 'inftimemonths')])
    known_times[, reps := .N, by = ptnum]
    known_times = known_times[!(is.na(inftimemonths) & reps > 1)]
    average_subset = merge(average_subset[, -c('incat_hiv', 'inftimemonths')], known_times[, -c('reps')], by = 'ptnum') 

    average_subset[, inftimeyears := inftimemonths/12]
}

# convert fragment number to numeric
average_subset[, fragment_int := as.numeric(substring(Fragment, 2))]

name = str_split(DATA_PATH, '/')[[1]][length(str_split(DATA_PATH, '/')[[1]])]
name = str_replace(name, '.csv', '')
dir.create(file.path(OUTPUT_PATH, 'adult_model_estimates'))
file_name = file.path(OUTPUT_PATH, 'adult_model_estimates', paste0('adult_model_estimates_', name, '.csv'))

fwrite(average_subset, file_name)

cat(file_name)
