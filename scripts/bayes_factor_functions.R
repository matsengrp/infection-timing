get_bayes_factor_results_path <- function(){
    path = file.path(OUTPUT_PATH, 'bayes_factor_results')
    dir.create(path, recursive = TRUE)
    filename = file.path(path, 'results.tsv')
    return(filename)
}

compile_results <- function(full_model, model_without_variation){
    filename = get_bayes_factor_results_path()

    # get marginal likelihood of original model
    full_model_log_marg_lik = bridge_sampler(full_model, maxiter = 200000, cores = NCPU, silent = TRUE)

    # get marginal likelihood of trial model
    model_without_var_log_marg_lik = bridge_sampler(model_without_variation, maxiter = 200000, cores = NCPU, silent = TRUE)

    # get bayes factor
    bayes_factor = bridgesampling::bf(full_model_log_marg_lik, model_without_var_log_marg_lik)

    results_df = data.table(original_model = MODEL_FILE, trial_model = NO_VAR_MODEL_FILE, original_model_log_marginal_likelihood_estimate = full_model_log_marg_lik[1], varying_variable = BAYES_FACTOR_VARIATION, original_model_bridge_sampling_iter = full_model_log_marg_lik[2], trial_model_log_marginal_likelihood_estimate = model_without_var_log_marg_lik[1], trial_model_bridge_sampling_iter = model_without_var_log_marg_lik[2], bayes_factor_result = bayes_factor[1]$bf)

    if (file.exists(filename)) {
        old_results = fread(filename)
        together = rbind(old_results, results_df)
        fwrite(together, filename, sep = '\t')
    } else {
        fwrite(results_df, filename, sep = '\t')
    }
    return(results_df)
}
