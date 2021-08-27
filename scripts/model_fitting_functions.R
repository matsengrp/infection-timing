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

check_data <- function(){
    raw_data = fread(INFANT_DATA_PATH)
    # check subject_ids
    raw_data[, extracted_subjects := as.numeric(search_sample_name_ptnum(Sample))]
    if (nrow(raw_data[extracted_subjects != ptnum]) > 0){
        print(paste0('Some subject ID mistakes found:'))
        print(raw_data[extracted_subjects != ptnum])
        print(paste0('Correcting subject ID mistakes now.'))
        raw_data[, ptnum := extracted_subjects]
    }
    # check fragment numbers
    raw_data[, extracted_frags := search_sample_name('F', Sample)]
    if (nrow(raw_data[extracted_frags != Fragment]) > 0){
        print(paste0('Some fragment ID mistakes found:'))
        print(raw_data[extracted_frags != Fragment])
        print(paste0('Correcting fragment ID mistakes now.'))
        raw_data[, Fragment := extracted_frags]
    }

    # get replicates
    raw_data[, replicate := search_sample_name('R', Sample)]
    return(raw_data[, -c('extracted_subjects', 'extracted_frags')])
}

preprocess_data <- function(data){
    # in cases where there are four replicates, filter to two
    data[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    data[, min_seq_run := min(seq_run), by = .(ptnum, Fragment, timepoint, replicate)]

    #TODO verify prioritization/filtering procedure
    # if there are more than 2 replicates, prioritize replicates from earlier runs 
    filtered = data[!(rep_count > 2 & min_seq_run != seq_run)]
    filtered[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    filtered = filtered[!(rep_count > 2 & replicate %in% c('R3', 'R4'))]
    filtered[, rep_count := .N, by = .(ptnum, Fragment, timepoint)]
    stopifnot(sum(unique(filtered$rep_count)) <= 3)
    return(filtered[, -c('rep_count', 'min_seq_run')])
}

index_subjects <- function(data){
    subjects = unique(data$ptnum)
    indices = seq(1:length(subjects))
    names(indices) = subjects
    data[, subject_index := indices[paste(ptnum)]]
    return(data)
}

configure_data <- function(){
    data = check_data()
    # remove NA cases
    #TODO verify NA removal
    data = data[!(is.na(pass2_APD))]
    if (isTRUE(PREPROCESS_DATA)){
        data = preprocess_data(data)
    }
    
    temp_cols = c('ptnum', 'month_visit', 'pass2_APD', 'Fragment', 'replicate', 'vload', 'incat_hiv')
    important_cols = temp_cols[!(temp_cols %in% c('replicate', 'pass2_APD'))]

    subset = data[, ..temp_cols]

    # average APD across replicates
    # TODO verify averaging procedure OR think about how to incorporate multiple measures in model?
    average_subset = subset[, mean(pass2_APD), by = important_cols]
    setnames(average_subset, 'V1', 'average_APD')

    # convert months to years
    average_subset[, year_visit := month_visit/12]
    # convert fragment number to numeric
    average_subset[, fragment_int := as.numeric(substring(Fragment, 2))]
    # assign index to each subject
    average_subset = index_subjects(average_subset)

    data_list = list(
                     observation_count = nrow(average_subset), 
                     subject_count = length(unique(average_subset$ptnum)), 
                     fragment_count = length(unique(average_subset$Fragment)), 
                     observed_time = average_subset$year_visit,
                     apd = average_subset$average_APD, 
                     subject = average_subset$subject_index,
                     subject_id = average_subset$ptnum,
                     fragment = average_subset$fragment_int,
                     infection_status = average_subset$incat_hiv
                     )
    return(data_list)
}

fit_model <- function(data, chains = 4, warmup_iterations = 10000, total_iterations = 20000, step_size = 0.99){
    model = stan(file = MODEL_FILE, 
                 data = data,
                 chains = chains, # number of Markov chains
                 warmup = warmup_iterations, # total number of warmup iterations per chain
                 iter = total_iterations, # total number of iterations per chain
                 cores = NCPU, # number of cores (can us up to one per chain)
                 control = list(adapt_delta = step_size, max_treedepth = 20)
                 )
    return(model)
}

