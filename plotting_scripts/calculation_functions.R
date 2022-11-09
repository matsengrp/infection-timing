calculate_mae_dt <- function(data, by_fragment = TRUE, var = 'mean_time_since_infection_estimate'){
    data[, abs_diff := abs(observed_time_since_infection - get(var))]
    if (isTRUE(by_fragment)){
        mae = data[, mean(abs_diff), by = fragment]
        setnames(mae, 'V1', 'mae')
        mae[, mae := round(mae, 3)]
    } else {
        mae = data[, mean(abs_diff)]
        mae = round(mae, 3)
    }
    return(mae)
}

calculate_rsme_dt <- function(data, frag){
    data[, abs_diff := abs(observed_time_since_infection - mean_time_since_infection_estimate)]
    data[, squared_diff := (observed_time_since_infection - mean_time_since_infection_estimate)^2]
    if (is.numeric(frag)){
        rmse = sqrt(sum(data[fragment == frag]$squared_diff)/length(data[fragment == frag]$squared_diff))
    } else {
        rmse = sqrt(sum(data$squared_diff)/length(data$squared_diff))
    }
    return(round(rmse, 3))
}

calculate_infection_time <- function(data_path, posterior_means){
    data = check_data(data_path, time_known = TRUE)
    # remove NA cases
    data = data[!(is.na(pass2_APD))]
    data = preprocess_data(data)
 
    post = unique(data[incat_hiv != 'IN UTERO'][, c('ptnum', 'inftimemonths')])
    post = post[!is.na(inftimemonths)]
    post[, it := inftimemonths/12]
    setnames(post, 'ptnum', 'subject_id')

    infant_data = configure_data(data_path)
    infant_data$number = seq(1, length(infant_data$is_post))
    infant_data = as.data.table(infant_data)

    subset = posterior_means[variable %like% 'correction']
    subset = subset[!(variable %like% 'utero')]
    subset$number = sapply(subset$variable, function(x) str_split(x, '\\[')[[1]][2])
    subset$number = as.numeric(str_remove(subset$number, '\\]'))
    subset = subset[order(number)]
    together = merge(subset, infant_data, by.x = 'number', by.y = 'subject')
    cols = c('5.5%', '94.5%', 'mean', 'infection_status', 'fragment', 'subject_id')
    subset2 = unique(together[, ..cols])
    setnames(subset2, 'mean', 'it')
    subset2 = subset2[infection_status == 'IN UTERO']
    subset2[, it := -1*it]
    return(unique(rbind(post[, c('subject_id', 'it')], subset2[, c('subject_id', 'it')])))
}
