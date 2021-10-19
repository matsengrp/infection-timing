data{
    int<lower=1> observation_count; // number of observations
    int<lower=1> subject_count;   // number of subjects
    int<lower=1> fragment_count;  // number of fragments
    real apd[observation_count];    // predictor variable
    int subject[observation_count];   // subject index
    int fragment[observation_count];  // fragment index
    int iteration_count;
    vector[iteration_count] subject_slope_delta_mean_estimate;
    vector[iteration_count] subject_slope_delta_variance_estimate;
    vector[iteration_count] baseline_slope;
    vector[iteration_count] fragment_slope_delta_reparameterized;
    vector[iteration_count] time_since_infection_variance_estimate;
} parameters {
} model {
}
generated quantities {
    matrix[iteration_count, observation_count] time_since_infection;
    matrix[iteration_count, subject_count] subject_slope_delta;
    matrix[iteration_count, observation_count] total_slope;
    matrix[iteration_count, observation_count] predicted_time_since_infection;
    for (subj in 1:subject_count){
        for (iter in 1:iteration_count){
            subject_slope_delta[iter, subj] = normal_rng( subject_slope_delta_mean_estimate[iter] , subject_slope_delta_variance_estimate[iter] );
          for (i in 1:observation_count){
              total_slope[iter, i] = (baseline_slope[iter] + subject_slope_delta[iter, subj] + fragment_slope_delta_reparameterized[iter]);
              predicted_time_since_infection[iter, i] = total_slope[iter, i] * apd[i]; // linear function relating apd to time_since_infection
              time_since_infection[iter, i] = normal_rng(predicted_time_since_infection[iter, i], time_since_infection_variance_estimate[iter]);
          }
        }
    }
}

