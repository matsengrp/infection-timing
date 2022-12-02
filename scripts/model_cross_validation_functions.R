reindex_subjects <- function(data){
    subjects = unique(data$subject_id)
    indices = seq(1:length(subjects))
    names(indices) = subjects
    data[, subject := indices[paste(subject_id)]]
    return(data)
}

reindex_fragments <- function(data){
    fragments = unique(data$fragment)
    indices = seq(1:length(fragments))
    names(indices) = fragments 
    data[, fragment_id := fragment]
    data[, fragment := indices[paste(fragment)]]
    return(data)
}

# configure already process data for use in leave one out 
configure_leave_one_out_data <- function(data_list, left_out_subject_id, model_posterior_for_prediction, type = "MCMC"){
    data = as.data.table(data_list) 

    # create left out subject data
    subset = data[subject_id == left_out_subject_id]
    subset = reindex_subjects(subset) 
    subset = reindex_fragments(subset)
    valid_fragments = unique(subset$fragment_id)
    if (type == 'MCMC') {
        temp_frag_param = model_posterior_for_prediction$fragment_slope_delta_reparameterized[, valid_fragments] 
    } else {
        temp_frag_param = model_posterior_for_prediction[names(model_posterior_for_prediction) %like% 'fragment_slope_delta_reparameterized']
        temp_frag_param = temp_frag_param[valid_fragments]
    }

    left_out_subject_data = list(
                     observation_count = nrow(subset), 
                     subject_count = length(unique(subset$subject_id)), 
                     fragment_count = length(unique(subset$fragment)), 
                     apd = subset$apd, 
                     subject = subset$subject,
                     subject_id = subset$subject_id,
                     fragment = subset$fragment,
                     fragment_id = subset$fragment_id,
                     iteration_count = length(model_posterior_for_prediction[['subject_slope_delta_mean_estimate']]),
                     subject_slope_delta_mean_estimate = model_posterior_for_prediction[['subject_slope_delta_mean_estimate']], 
                     subject_slope_delta_variance_estimate = model_posterior_for_prediction[['subject_slope_delta_variance_estimate']],
                     baseline_slope = model_posterior_for_prediction[['baseline_slope']],
                     fragment_slope_delta_reparameterized = temp_frag_param,
                     time_since_infection_variance_estimate = model_posterior_for_prediction[['time_since_infection_variance_estimate']]
                     )
    return(left_out_subject_data)
}

configure_data_without_left_out_subject <- function(data_list, left_out_subject_id){
    data = as.data.table(data_list)
    # create data for all other subjects 
    other_data = data[subject_id != left_out_subject_id]
    other_data = reindex_subjects(other_data)
    other_subject_data = list(
        observation_count = nrow(other_data), 
        subject_count = length(unique(other_data$subject_id)), 
        fragment_count = length(unique(other_data$fragment)), 
        observed_time_since_infection = other_data$observed_time_since_infection,
        apd = other_data$apd, 
        subject = other_data$subject,
        subject_id = other_data$subject,
        fragment = other_data$fragment,
        infection_status = other_data$infection_status, 
        is_utero = other_data$is_utero,
        is_post = other_data$is_post
    )
    return(other_subject_data) 
}

run_leave_one_out_validation <- function(data, type = 'MCMC'){
    stopifnot(type %in% c('MCMC', 'MAP'))
    if (type == 'MCMC') {
        results = run_leave_one_out_validation_mcmc(data)
    } else if (type == 'MAP') {
        results = run_leave_one_out_validation_map(data)
    } 
    return(results)
}

run_leave_one_out_validation_map <- function(data){
    together = data.table() 
    for (subject in unique(data$subject_id)){
        # remove left out subject from data
        other_subjects = configure_data_without_left_out_subject(data, subject)
        # fit left out model
        other_subject_model = fit_model(other_subjects, type = 'MAP')
        # get posteriors
        model_posterior_for_prediction = other_subject_model$par
        # get data for left out subject
        leave_one_out_data = configure_leave_one_out_data(data, subject, model_posterior_for_prediction, type = 'MAP')
        # get predictions
        leave_one_out_predictions = run_prediction(leave_one_out_data, type = 'MAP')
        # get prediction posteriors
        leave_one_out_posteriors = get_prediction_posterior(leave_one_out_predictions, type = 'MAP')

        # clean data
        print(paste0('finished validation for subject ', subject))

        transformed_posteriors = transform_posterior_matrix_to_dataframe(leave_one_out_data, leave_one_out_posteriors)
        together = rbind(together, transformed_posteriors)
    }
    # stopImplicitCluster()
    together$time_correction_type = TIME_CORRECTION_TYPE
    return(together)
}

run_leave_one_out_validation_mcmc <- function(data){
    together = data.table() 
    # registerDoParallel(cores=NCPU)
    # together = foreach(subject = unique(data$subject_id), .combine = 'rbind') %dopar% {
    for (subject in unique(data$subject_id)){
        # remove left out subject from data
        other_subjects = configure_data_without_left_out_subject(data, subject)
        # fit left out model
        other_subject_model = fit_model(other_subjects, chains = 1)
        # get posteriors
        model_posterior_for_prediction = rstan::extract(other_subject_model)
        # get data for left out subject
        leave_one_out_data = configure_leave_one_out_data(data, subject, model_posterior_for_prediction, type = 'MCMC')
        # get predictions
        leave_one_out_predictions = run_prediction(leave_one_out_data)
        # get prediction posteriors
        leave_one_out_posteriors = get_prediction_posterior(leave_one_out_predictions)
        # clean data
        print(paste0('finished validation for subject ', subject))

        transformed_posteriors = transform_posterior_matrix_to_dataframe(leave_one_out_data, leave_one_out_posteriors)
        together = rbind(together, transformed_posteriors)
    }
    # stopImplicitCluster()
    together$time_correction_type = TIME_CORRECTION_TYPE
    return(together)
}
