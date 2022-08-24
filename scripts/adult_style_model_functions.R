get_adult_style_model_name <- function(fragment){
    path = file.path(PROJECT_PATH, 'scripts', 'adult_style_models', 'model_fits')
    dir.create(path, recursive = TRUE)
    data_description = str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]][length(str_split(TRAINING_INFANT_DATA_PATH, '/')[[1]])]
    data_description = str_split(data_description, '.csv')[[1]][1]
    model_name = paste0('/fragment_', fragment, '_adult_style_model_', data_description, '.rds')
    together = paste0(path, model_name)
    return(together)
}

get_adult_style_model_predictions_name <- function(data_name){
    output_path = file.path(OUTPUT_PATH, 'adult_style_model_predictions')
    dir.create(output_path, recursive = TRUE)
    data_name = str_remove(data_name, '.csv')
    data_name = str_remove(data_name, '.tsv')
    data_name = str_split(data_name, '/')[[1]][length(str_split(data_name, '/')[[1]])]
    file_name = paste0(data_name, '_adult_style_model_results.tsv')
    together = file.path(output_path, file_name)
    return(together)
}

get_adult_style_model_loocv_name <- function(){
    output_path = file.path(OUTPUT_PATH, 'adult_style_model_loocv')
    dir.create(output_path, recursive = TRUE)
    data_name = str_remove(TRAINING_INFANT_DATA_PATH, '.csv')
    data_name = str_remove(data_name, '.tsv')
    data_name = str_split(data_name, '/')[[1]][length(str_split(data_name, '/')[[1]])]
    file_name = paste0('loocv_', data_name, '_adult_style_model_results.tsv')
    together = file.path(output_path, file_name)
    return(together)
}

run_adult_style_model_loocv <- function(processed_training_data){
    processed_training_data_dt = as.data.table(data.frame(processed_training_data))

    results = data.table()
    for (frag in unique(processed_training_data_dt$fragment)){
        fragment_processed_training_data_dt = processed_training_data_dt[fragment == frag]
        for (infant in unique(fragment_processed_training_data_dt$subject)){
            temp_data = fragment_processed_training_data_dt[subject != infant]
            held_out_data = fragment_processed_training_data_dt[subject == infant]
            if (nrow(held_out_data) > 0){
                temp_model = lad(observed_time_since_infection ~ apd, data = temp_data, na.action = na.pass)
                held_out_data$adult_style_model_time_since_infection_estimates = stats::predict(temp_model, newdata = held_out_data, na.action = na.pass)
                held_out_data$model_type = paste0(frag, '_fragment_model')
                results = rbind(results, held_out_data)
            }
        }
    }
    return(results)
}
