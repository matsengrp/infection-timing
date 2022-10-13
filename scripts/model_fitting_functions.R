source(paste0(PROJECT_PATH, '/scripts/model_time_correction_type/', TIME_CORRECTION_TYPE, '.R'))

search_sample_name <- function(char, column){
    split_col = str_split(column, '_')
    char_list = unlist(purrr::map_depth(split_col, 1, function(x) grep(char, x, value = TRUE)))
    return(char_list)
}

search_sample_name_ptnum <- function(column){
    rm_pt = str_remove_all(column, 'pt')
    split_col = str_split(rm_pt, '_')
    char_list = unlist(purrr::map_depth(split_col, 1, function(x) sub(".*?(\\d+).*", "\\1", paste(x, collapse = ' ')))) 
    return(char_list)
}

check_data <- function(data_path, time_known = TRUE){
    raw_data = fread(data_path)
    if (isTRUE(time_known)) {
        temp_cols = c('Sample', 'ptnum', 'month_visit', 'pass2_APD', 'Fragment', 'vload', 'incat_hiv')
    } else {
        temp_cols = c('Sample', 'ptnum', 'pass2_APD', 'Fragment')
    }
    stopifnot(all(temp_cols %in% colnames(raw_data)))

    # check subject_ids
    raw_data[!(Sample == ''), extracted_subjects := as.numeric(search_sample_name_ptnum(Sample))]
    raw_data[Sample == '', extracted_subjects := ptnum]
    if (nrow(raw_data[extracted_subjects != ptnum]) > 0){
        raw_data[, ptnum := extracted_subjects]
    }
    # check fragment numbers
    raw_data[!(Sample == ''), extracted_frags := search_sample_name('F', Sample)]
    raw_data[Sample == '', extracted_frags := Fragment]
    if (nrow(raw_data[extracted_frags != Fragment]) > 0){
        raw_data[, Fragment := extracted_frags]
    }

    # filter out no APD cases
    raw_data = raw_data[!is.na(pass2_APD)]
    # get replicates
    raw_data[, replicate := search_sample_name('R', Sample)]
    return(raw_data[, -c('extracted_subjects', 'extracted_frags')])
}

preprocess_data <- function(data){
    # in cases where there are four replicates, filter to two
    data[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    data[, min_seq_run := min(seq_run), by = .(ptnum, Fragment, timepoint, replicate)]
    
    # fix duplicates
    data = unique(data[, -c('vload')])
    
    # if there are more than 2 replicates, prioritize replicates from earlier runs 
    filtered = data[!(rep_count > 2 & min_seq_run != seq_run)]
    filtered[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    filtered = filtered[!(rep_count > 2 & replicate %in% c('R3', 'R4'))]
    filtered[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    filtered = filtered[!(rep_count > 2 & replicate %in% c('R5', 'R6'))]
    filtered[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    stopifnot(sum(unique(filtered$rep_count)) <= 3)
    return(filtered[, -c('rep_count', 'min_seq_run')])
}

filter_data <- function(data){
    # filter subjects/fragment pairs for which we only have one sample timepoint
    data[, timepoint_count := .N, by = .(ptnum, Fragment)]
    filtered = data[timepoint_count > 1]
    return(filtered[, -c('timepoint_count')])
}

index_subjects <- function(data){
    subjects = unique(data$ptnum)
    indices = seq(1:length(subjects))
    names(indices) = subjects
    data[, subject_index := indices[paste(ptnum)]]
    return(data)
}

index_infection_time <- function(data){
    it = unique(data$incat_hiv)
    indices = seq(1:length(it))
    names(indices) = it
    data[incat_hiv == it[1], infection_time_index := indices[1]]
    data[incat_hiv == it[2], infection_time_index := indices[2]]
    stopifnot(length(it) == 2)
    return(data)
}


configure_data <- function(data_path, time_known = TRUE){
    data = check_data(data_path, time_known)
    # remove NA cases
    data = data[!(is.na(pass2_APD))]
    if (isTRUE(PREPROCESS_DATA)){
        data = preprocess_data(data)
    }

    temp_cols = c('ptnum', 'month_visit', 'pass2_APD', 'Fragment', 'replicate', 'incat_hiv', 'inftimemonths')

    important_cols = temp_cols[!(temp_cols %in% c('replicate', 'pass2_APD'))]

    subset = data[, ..temp_cols]

    # average APD across replicates
    average_subset = subset[, mean(pass2_APD), by = important_cols]
    setnames(average_subset, 'V1', 'average_APD')

    # fill in missing infection times
    known_times = unique(subset[, c('ptnum', 'incat_hiv', 'inftimemonths')])
    known_times[, reps := .N, by = ptnum]
    known_times = known_times[!(is.na(inftimemonths) & reps > 1)]
    average_subset = merge(average_subset[, -c('incat_hiv', 'inftimemonths')], known_times[, -c('reps')]) 

    # filter out subjects, fragments without multiple timepoints
    # average_subset = filter_data(average_subset)

    # convert months to years
    average_subset[, year_visit := month_visit/12]
    average_subset[, inftimeyears := inftimemonths/12]
    # convert fragment number to numeric
    average_subset[, fragment_int := as.numeric(substring(Fragment, 2))]
    # assign index to each subject
    average_subset = index_subjects(average_subset)
    average_subset = index_infection_time(average_subset)
    # for late infected subjects, convert observed time to time since infection given estimated infection time
    average_subset[incat_hiv != 'IN UTERO', year_visit := year_visit-inftimeyears]

    average_subset[incat_hiv == 'IN UTERO', is_utero := TRUE]
    average_subset[incat_hiv != 'IN UTERO', is_utero := FALSE]
    average_subset[incat_hiv == 'LATE, after M1', is_post := TRUE]
    average_subset[incat_hiv != 'LATE, after M1', is_post := FALSE]

    data_list = list(
                     observation_count = nrow(average_subset), 
                     subject_count = length(unique(average_subset$ptnum)), 
                     fragment_count = length(unique(average_subset$Fragment)), 
                     infection_time_count = length(unique(average_subset$incat_hiv)),
                     observed_time_since_infection = average_subset$year_visit,
                     apd = average_subset$average_APD, 
                     subject = average_subset$subject_index,
                     subject_id = average_subset$ptnum,
                     fragment = average_subset$fragment_int,
                     infection_time = average_subset$infection_time_index, 
                     infection_status = average_subset$incat_hiv, 
                     is_utero = average_subset$is_utero,
                     is_post = average_subset$is_post
                     )
    return(data_list)
}

fit_model <- function(data, chains = 4, warmup_iterations = 10000, total_iterations = 20000, step_size = 0.99, model_file = MODEL_FILE, type = 'MCMC'){
    stopifnot(type %in% c('MCMC', 'MAP'))
    if (type == 'MCMC') {
        model = fit_model_mcmc(data, chains, warmup_iterations, total_iterations, step_size, model_file)
    } else if (type == 'MAP') {
        model_code = stan_model(model_file)
        model = optimizing(model_code, data, hessian = TRUE)
    } 
    return(model)
}

fit_model_mcmc <- function(data, chains = 4, warmup_iterations = 10000, total_iterations = 20000, step_size = 0.99, model_file = MODEL_FILE){
    model = stan(file = model_file, 
                 data = data,
                 chains = chains, # number of Markov chains
                 warmup = warmup_iterations, # total number of warmup iterations per chain
                 iter = total_iterations, # total number of iterations per chain
                 cores = NCPU, # number of cores (can us up to one per chain)
                 seed = 555,
                 control = list(adapt_delta = step_size, max_treedepth = 20)
                 )
    return(model)
}

get_model_fit_name <- function(type){
    path = file.path(PROJECT_PATH, 'scripts', 'stan_models', 'model_fits')
    dir.create(file.path(PROJECT_PATH, 'scripts', 'stan_models', 'cred_intervals'), recursive = TRUE)
    dir.create(path, recursive = TRUE)
    name = str_split(MODEL_FILE, '/')[[1]][7]
    name = str_split(name, '.stan')[[1]][1]
    data_description = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
    data_description = str_split(data_description, '.csv')[[1]][1]
    model_name = paste0('/', name, '_', data_description, '_', type, '.rds')
    together = paste0(path, model_name)
    return(together)
}

save_model_fit <- function(model, model_name = get_model_fit_name('MCMC')){
    saveRDS(model, model_name)
}

load_model_fit <- function(model_name = get_model_fit_name('MCMC')){
    model = readRDS(model_name)
    return(model)
}


