library(rethinking)
library(data.table)
library(dplyr)

infant_data = read.csv("../_ignore/AllRunsAvg.csv")

preprocess_noavg <- function(data){
    data = data.table(data)
    data_avg = data[APD_1 != 'NA']
    sample_name_list = strsplit(as.character(data_avg$SampleName), "_")
    sample_name_list_updated = c()
    for (i in seq(length(sample_name_list))){
        if (sample_name_list[[i]][1] %in% c("pt93", "pt108", "pt135", "pt170", "pt45")){
            sample_name_list_updated = c(sample_name_list_updated, sample_name_list[[i]][3])
        } else if (sample_name_list[[i]][1] %in% c("pt258", "pt313", "pt485", "pt490")){
            sample_name_list_updated = c(sample_name_list_updated, sample_name_list[[i]][4])
        } else if (sample_name_list[[i]][1] %in% c("pt", "Pt")){
            sample_name_list_updated = c(sample_name_list_updated, sample_name_list[[i]][5])
        }
    }     
    data_avg2 = data_avg %>% select(Sample, ActualTOI..year., Fragment, VL, APD_1, APD_5, APD_10)
    colnames(data_avg2) = c("subject", "time", "fragment", "vload","avg_apd1", "avg_apd5", "avg_apd10")
    data_avg2$run = sample_name_list_updated
    data_avg2$run[data_avg2$run == 'R2'] <- 2
    data_avg2$run[data_avg2$run == 'R1'] <- 1
#    data_avg3 = data_avg2[, multiple_run_index := .GRP, by = .(run, fragment)]
#    data_avg2$time = data_avg2$time + 0.125
    return(data_avg2)
}

save_models_plots_stats <- function(model, type){
    if (type == "fragment_subject_varying_slopes"){
        intercepts = c('total_intercept')
        slopes = c('subject_slope','fragment_slope','population_avg_slope')
    } else if (type == "fragment_subject_varying_slopes_no_int"){
        slopes = c('subject_slope','fragment_slope','population_avg_slope')
    } else if (type == "time_error"){
        slopes = c('subject_slope','fragment_slope','population_avg_slope')
    } else if (type == "time_error_int"){
        intercepts = c('total_intercept')
        slopes = c('subject_slope','fragment_slope','population_avg_slope')
    } else {
        print('ERROR')
    }

    if (type == "subject_varying_intercepts_slopes" | type == "subject_varying_intercepts_fragment_subject_varying_slopes"| type == "subject_varying_intercepts_fragment_subject_run_varying_slopes" | type == "time_error"| type == "time_error_int"){
        print("no pairs plot")
    } else {
        pdf(paste0("plots/",model,"/pairs_plot.pdf"))
        pairs(get(model))
        dev.off()
    }

    if (length(slopes) == 1){
        print('only one slope')
    } else{
        pdf(paste0("plots/",model,"/slope_plots.pdf"))
        plot(precis(get(model), depth = 2, pars= paste(slopes)))
        dev.off()
    }
    
    # plot trace plot and posterior distribution to make sure the effective number of subjects and Rhat values look alright
    pdf(paste0("plots/",model,"/convergence_plots.pdf"))
    plot(get(model))
    dev.off()

    # Save model
    saveRDS(get(model), paste0("models/", model, ".rds"))

    # Look at model statistics
    print(precis(get(model), depth = 2))
}

save_models <- function(model){
    # Save model
    saveRDS(get(model), paste0("models/", model, ".rds"))
}

## Clean data and index subjects and fragments
index_data_subjects_frags <- function(data, apd){
    data_cleaned = preprocess_noavg(data)
    data_cleaned$apd = data_cleaned[[paste0('avg_apd', apd)]]
    data_cleaned = data_cleaned[order(subject, fragment)]
    subject_index = c(seq(1, length(unique(data_cleaned$subject))))
    names(subject_index) = unique(data_cleaned$subject)
    frag_index = c(seq(1, length(unique(data_cleaned$fragment))))
    names(frag_index) = unique(data_cleaned$fragment)
    indexed_data_cleaned = data.frame()
    for (samp in unique(data_cleaned$subject)){
        data_cleaned1 = NULL
        data_cleaned1 = data_cleaned[with(data_cleaned, subject == samp)]
        data_cleaned1$subject_index = subject_index[[samp]]
        for (frag in unique(data_cleaned$fragment)){
            data_cleaned2 = NULL
            data_cleaned2 = data_cleaned1[with(data_cleaned1, fragment == frag)]
            data_cleaned2$frag_index = frag_index[[frag]]
            indexed_data_cleaned = rbind(indexed_data_cleaned, data_cleaned2)
        }   
    }
    return(indexed_data_cleaned)
}

# indexed infant data
indexed_infant_data_cleaned = index_data_subjects_frags(infant_data, apd = 1)

# Make data.table with only necessary columns
necessary_columns <- c("time", "apd", "subject_index", "frag_index")
indexed_infant_data_cleaned_subset = indexed_infant_data_cleaned[, ..necessary_columns]

indexed_infant_data_cleaned_subject_subset = index_data_subjects_frags(data.table(infant_data)[!(Sample %in% c('pt485', 'pt258', 'pt313'))], apd = 1)

indexed_infant_data_cleaned_subset_subject_subset = indexed_infant_data_cleaned_subject_subset[, ..necessary_columns]