configure_newdata <- function(data_path, time_known = FALSE){
    data = check_data(data_path, time_known)
    data = data[!(is.na(pass2_APD))]

    temp_cols = c('ptnum', 'pass2_APD', 'Fragment')
    important_cols = temp_cols[!(temp_cols %in% c('replicate', 'pass2_APD'))]

    subset = data[, ..temp_cols]

    # average APD across replicates
    average_subset = subset[, mean(pass2_APD), by = important_cols]
    setnames(average_subset, 'V1', 'average_APD')

    # convert fragment number to numeric
    average_subset[, fragment_int := as.numeric(substring(Fragment, 2))]
    # assign index to each subject
    average_subset = index_subjects(average_subset)

    data_list = list(
                     observation_count = nrow(average_subset), 
                     subject_count = length(unique(average_subset$ptnum)), 
                     fragment_count = length(unique(average_subset$Fragment)), 
                     apd = average_subset$average_APD, 
                     subject = average_subset$subject_index,
                     subject_id = average_subset$ptnum,
                     fragment = average_subset$fragment_int
                     )
    return(data_list)
}

# configure already process data for use in leave one out 
configure_data_for_prediction <- function(data_list, model_posterior_for_prediction, type = 'MCMC'){
    data = as.data.table(data_list) 
    data$fragment_id = data$fragment
    valid_fragments = unique(data$fragment_id)
    if (type == 'MCMC') {
        temp_frag_param = model_posterior_for_prediction$fragment_slope_delta_reparameterized[, valid_fragments] 
    } else {
        temp_frag_param = model_posterior_for_prediction[names(model_posterior_for_prediction) %like% 'fragment_slope_delta_reparameterized']
    }
    subject_data = list(
                     observation_count = nrow(data), 
                     subject_count = length(unique(data$subject_id)), 
                     fragment_count = length(unique(data$fragment)), 
                     apd = data$apd, 
                     subject = data$subject,
                     subject_id = data$subject_id,
                     fragment = data$fragment,
                     fragment_id = data$fragment_id,
                     iteration_count = length(model_posterior_for_prediction[['subject_slope_delta_mean_estimate']]),
                     subject_slope_delta_mean_estimate = model_posterior_for_prediction[['subject_slope_delta_mean_estimate']], 
                     subject_slope_delta_variance_estimate = model_posterior_for_prediction[['subject_slope_delta_variance_estimate']],
                     baseline_slope = model_posterior_for_prediction[['baseline_slope']],
                     fragment_slope_delta_reparameterized = temp_frag_param,
                     time_since_infection_variance_estimate = model_posterior_for_prediction[['time_since_infection_variance_estimate']]
                     )
    return(subject_data)
}

get_prediction_stan_file <- function(type, frag_count){
    base_file = paste0(PROJECT_PATH, '/scripts/stan_models/')
    name = paste0('predict.stan')
    if (frag_count == 1){
        name = str_replace(name, 'predict', 'predict_one_frag')
    }
    if (type != 'MCMC'){
        name = str_replace(name, '.stan', '_map.stan')
    }
    tog = paste0(base_file, name)
    return(tog)
}

run_prediction <- function(data, total_iterations = 100, type = 'MCMC'){
    frag_count = data$fragment_count
    prediction_file = get_prediction_stan_file(type, frag_count)
    predictions = stan(file = prediction_file, 
                       data = data,
                       iter = total_iterations, # total number of iterations per chain
                       seed = 555,
                       algorithm = 'Fixed_param'
                       )
    return(predictions)
}

get_prediction_posterior <- function(predictions, type = 'MCMC'){
    prediction_posteriors = rstan::extract(predictions)
    if (type == 'MCMC') {
        iterations = dim(prediction_posteriors$time_since_infection)[2]
        observations = dim(prediction_posteriors$time_since_infection)[3]
        time_posteriors = matrix(NA, nrow = iterations, ncol = observations)
        for (iter in 1:iterations) {
            for (obs in 1:observations) {
                time_posteriors[iter, obs] = median(prediction_posteriors$time_since_infection[, iter, obs])
            }
        }
    } else if (type == 'MAP'){
        iterations = 1
        observations = dim(prediction_posteriors$time_since_infection)[2]
        time_posteriors = matrix(NA, nrow = iterations, ncol = observations)
        for (iter in 1:iterations) {
            for (obs in 1:observations) {
                time_posteriors[iter, obs] = median(prediction_posteriors$time_since_infection[,obs])
            }
        }
    }
    return(time_posteriors)
}

transform_posterior_matrix_to_dataframe <- function(data, posteriors) {
    dt = as.data.table(data)
    data_cols = c('subject_id', 'fragment_id', 'apd')
    dt = unique(dt[, ..data_cols])
    dt$apd = as.character(dt$apd)

    posterior_dt = as.data.table(posteriors)
    colnames(posterior_dt) = dt$apd
    posterior_dt$iteration = seq(1:nrow(posterior_dt))

    posterior_dt = posterior_dt %>% pivot_longer(!iteration, names_to = 'apd', values_to = 'time_since_infection_posterior_draw') %>% as.data.table()

    together = merge(dt, posterior_dt, by = 'apd')
    together$apd = as.numeric(together$apd)
    return(together)
}

predict_posterior <- function(data, model, newdata = TRUE, type = 'MCMC'){
    if (type == 'MCMC') {
        posterior = rstan::extract(model)
    } else {
        posterior = model$par
    }
    if (isTRUE(newdata)){
        configured_data = configure_data_for_prediction(data, posterior, type)
        # get predictions
        predictions = run_prediction(configured_data, type)
        # get prediction posteriors
        prediction_posteriors = get_prediction_posterior(predictions, type)
    } else {
        prediction_posteriors = posterior$predicted_time_since_infection    
    }
    transformed_posteriors = transform_posterior_matrix_to_dataframe(data, prediction_posteriors)
    medians = get_posterior_means(transformed_posteriors, data)
    return(list(posteriors = transformed_posteriors, medians = medians))
}

get_mean_observed_time_to_time_since_infection_prediction <- function(data, model, type = "MCMC"){
    if (type == 'MCMC') {
        posterior = rstan::extract(model)
    } else {
        posterior = model$par
    }
    time_correction_posteriors = data.table(posterior[['observed_time_to_time_since_infection_correction']])
    colnames(time_correction_posteriors) = unique(paste0('subject_', data$subject_id))

    means = time_correction_posteriors[, lapply(.SD, mean)]
    means = means %>% 
        pivot_longer(everything(), names_to = 'subject_id', values_to = 'observed_time_correction') %>%
        as.data.table()

    means$subject_id = substring(means$subject_id, 9)
    return(means)
}

get_posterior_means <- function(posteriors, all_subject_data){
    posterior_means = posteriors[, median(time_since_infection_posterior_draw), by = .(subject_id, fragment_id, apd)]
    setnames(posterior_means, 'V1', 'mean_time_since_infection_estimate')
    setnames(posterior_means, 'fragment_id', 'fragment')
    all_subject_data = as.data.table(all_subject_data)
    all_subject_data$apd = as.numeric(as.character(all_subject_data$apd))
    merged = merge(all_subject_data, posterior_means, by = c('subject_id', 'fragment', 'apd'))

    stopifnot(nrow(merged) == nrow(all_subject_data))
    return(merged)
}

predict <- function(data, model, type = 'MCMC'){
    newdata = ifelse(TRAINING_INFANT_DATA_PATH == TESTING_INFANT_DATA_PATH, FALSE, TRUE)

    posteriors = predict_posterior(data, model, newdata, type)
    posteriors = posteriors$posteriors
    posterior_mean_dt = posteriors$medians
    #TODO alter this once using complete dataset
    if (isFALSE(newdata)){
        observed_time_corrections = get_mean_observed_time_to_time_since_infection_prediction(data, model, type)
        posterior_mean_dt = merge(posterior_mean_dt, observed_time_corrections, by = 'subject_id')
        posterior_mean_dt$corrected_observed_time = as.numeric(posterior_mean_dt$observed_time) + as.numeric(posterior_mean_dt$observed_time_correction)
    }

    return(posterior_mean_dt)
}

simulate_apd_time_stan <- function(model){
    apd = seq(0.0001, 0.025, by = 0.0004) 
    post = rstan::extract(model)

    simulations = foreach(index = seq(30), .combine = 'rbind') %do% {
        simulate_subject = rnorm(1,post$subject_slope_delta_mean_estimate[index], post$subject_slope_delta_variance_estimate[index])
        simulate_fragment = rnorm(1,post$fragment_slope_delta_mean_estimate[index], post$fragment_slope_delta_variance_estimate[index])
        baseline_slope = post$baseline_slope[index]
        time_variance = post$time_since_infection_variance_estimate[index]
        time_since_infection = rnorm(length(apd), (simulate_subject + simulate_fragment + baseline_slope) * apd, time_variance)
        data.table(apd = apd, predicted_time_since_infection = time_since_infection)
    }
    return(simulations)
}

get_model_output_filename <- function(input_data, type = 'MCMC'){
    path = file.path(OUTPUT_PATH, 'model_predictions')
    if (type != 'MCMC'){
        path = file.path(path, type)
    }
    dir.create(path, recursive = TRUE)
    data_name = str_remove(input_data, '.csv')
    data_name = str_remove(data_name, '.tsv')
    data_name = str_split(data_name, '/')[[1]][length(str_split(data_name, '/')[[1]])]
    model_name = str_split(MODEL_FILE, '/')[[1]][length(str_split(MODEL_FILE, '/')[[1]])]
    model_name = str_remove(model_name, '.stan')
    file_name = paste0(data_name, '_', model_name, '_results.tsv')
    together = file.path(path, file_name)
    return(together)
}

save_prediction_results <- function(results, input_data, type = 'MCMC', file_name = get_model_output_filename(input_data, type)){ 
    fwrite(results, file_name, sep = '\t')
}

get_posterior_interval_data <- function(model){
    matrix_data = as.matrix(model)
    require(rstanarm)
    post = posterior_interval(matrix_data)
    post_df = as.data.frame(post)
    post_df$variable = rownames(post_df)
    post_dt = as.data.table(post_df)
    mean = apply(matrix_data, 2, mean)
    mean_df = as.data.frame(mean)
    mean_df$variable = rownames(mean_df)
    mean_dt = as.data.table(mean_df)
    tog = merge(post_dt, mean_dt, by = 'variable')
    return(tog)
}

get_posterior_prediction_coverage <- function(loocv_results, actual_times){
    quants = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
    for (q in quants){
        loocv_results[, paste0('q', q) := quantile(time_since_infection_posterior_draw, probs = q), by = .(subject_id, fragment_id, apd)]
    }
    cols = c('apd', 'subject_id', 'fragment_id', 'q0', 'q0.05', 'q0.25', 'q0.5', 'q0.75', 'q0.95', 'q1')
    results = unique(loocv_results[, ..cols])
    setnames(results, 'fragment_id', 'fragment')
    results$apd = round(results$apd, 8)
    actual_times$apd = round(actual_times$apd, 8)
    tog = merge(results, actual_times, by = c('subject_id', 'apd', 'fragment'))
    # collapse all fragments
    quant_cols = paste0('q', quants)
    cols = c('subject_id', quant_cols, 'observed_time_since_infection')
    subset = unique(tog[, ..cols])
    subset = subset[,lapply(.SD, mean), by=.(subject_id,observed_time_since_infection)]
    # calculate coverage
    subset[observed_time_since_infection > q0.05 & observed_time_since_infection < q0.95, coverage_90 := TRUE]
    subset[!(observed_time_since_infection > q0.05 & observed_time_since_infection < q0.95), coverage_90 := FALSE]
    subset[observed_time_since_infection > q0.25 & observed_time_since_infection < q0.75, coverage_50 := TRUE]
    subset[!(observed_time_since_infection > q0.25 & observed_time_since_infection < q0.75), coverage_50 := FALSE]
    coverage_50 = subset[, .N, by = coverage_50][, prop := N/sum(N)]
    coverage_90 = subset[, .N, by = coverage_90][, prop := N/sum(N)]
    return(list(data = subset, coverage_50 = coverage_50, coverage_90 = coverage_90))
}
