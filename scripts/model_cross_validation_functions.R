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
configure_leave_one_out_data <- function(data_list, left_out_subject_id, model_posterior_for_prediction){
    data = as.data.table(data_list) 

    # create left out subject data
    subset = data[subject_id == left_out_subject_id]
    subset = reindex_subjects(subset) 
    subset = reindex_fragments(subset)
    valid_fragments = unique(subset$fragment_id)
    left_out_subject_data = list(
                     observation_count = nrow(subset), 
                     subject_count = length(unique(subset$subject_id)), 
                     fragment_count = length(unique(subset$fragment)), 
                     apd = subset$apd, 
                     subject = subset$subject,
                     subject_id = subset$subject_id,
                     fragment = subset$fragment,
                     fragment_id = subset$fragment_id,
                     iteration_count = length(model_posterior_for_prediction$subject_slope_delta_mean_estimate),
                     subject_slope_delta_mean_estimate = model_posterior_for_prediction$subject_slope_delta_mean_estimate, 
                     subject_slope_delta_variance_estimate = model_posterior_for_prediction$subject_slope_delta_variance_estimate,
                     baseline_slope = model_posterior_for_prediction$baseline_slope,
                     fragment_slope_delta_reparameterized = model_posterior_for_prediction$fragment_slope_delta_reparameterized[,valid_fragments],
                     time_since_infection_variance_estimate = model_posterior_for_prediction$time_since_infection_variance_estimate
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

run_prediction <- function(data, chains = 1, total_iterations = 100){
    prediction_file = paste0(PROJECT_PATH, '/scripts/stan_models/predict.stan')
    predictions = stan(file = prediction_file, 
                 data = data,
                 chains = chains, # number of Markov chains
                 iter = total_iterations, # total number of iterations per chain
                 cores = NCPU, # number of cores (can us up to one per chain)
                 seed = 555,
                 algorithm = 'Fixed_param'
                 )
    return(predictions)
}

get_leave_one_out_prediction_posterior <- function(leave_one_out_data, predictions){
    prediction_posteriors = rstan::extract(predictions)
    iterations = dim(prediction_posteriors$time_since_infection)[2]
    observations = dim(prediction_posteriors$time_since_infection)[3]
    time_posteriors = matrix(NA, nrow = iterations, ncol = observations)
    for (iter in 1:iterations) {
        for (obs in 1:observations) {
            time_posteriors[iter, obs] = mean(prediction_posteriors$time_since_infection[, iter, obs])
        }
    }
    return(time_posteriors)
}

transform_posterior_matrix_to_dataframe <- function(leave_one_out_data, posteriors) {
    dt = as.data.table(leave_one_out_data)
    data_cols = c('subject_id', 'fragment_id', 'apd')
    dt = unique(dt[, ..data_cols])
    dt$apd = as.character(dt$apd)

    posterior_dt = as.data.table(posteriors)
    colnames(posterior_dt) = dt$apd
    posterior_dt$iteration = seq(1:nrow(posterior_dt))

    posterior_dt = posterior_dt %>% pivot_longer(!iteration, names_to = 'apd', values_to = 'time_since_infection_posterior_draw') %>% as.data.table()

    together = merge(dt, posterior_dt, by = 'apd')
    together$apd = as.numeric(together$ap)
    return(together)
}

run_leave_one_out_validation <- function(data){
    together = data.table() 
    for (subject in unique(data$subject_id)){
        # remove left out subject from data
        other_subjects = configure_data_without_left_out_subject(data, subject)
        # fit left out model
        other_subject_model = fit_model(other_subjects, chains = 1)
        # get posteriors
        model_posterior_for_prediction = rstan::extract(other_subject_model)
        # get data for left out subject
        leave_one_out_data = configure_leave_one_out_data(data, subject, model_posterior_for_prediction)
        # get predictions
        leave_one_out_predictions = run_prediction(leave_one_out_data)
        # get prediction posteriors
        leave_one_out_posteriors = get_leave_one_out_prediction_posterior(leave_one_out_data, leave_one_out_predictions)
        # clean data
        leave_one_out_posteriors_dt = transform_posterior_matrix_to_dataframe(leave_one_out_data, leave_one_out_posteriors)
        together = rbind(together, leave_one_out_posteriors_dt)
    }
    together$time_correction_type = TIME_CORRECTION_TYPE
    return(together)
}
