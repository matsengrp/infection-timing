calculate_mae_dt <- function(data){
    data[, abs_diff := abs(observed_time_since_infection - mean_time_since_infection_estimate)]
    mae = data[, mean(abs_diff), by = fragment]
    setnames(mae, 'V1', 'mae')
    mae[, mae := round(mae, 3)]
    return(mae)
}
