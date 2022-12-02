data{
    int<lower=1> observation_count; // number of observations
    int<lower=1> subject_count;   // number of subjects
    int<lower=1> fragment_count;  // number of fragments
    real apd[observation_count];    // predictor variable
    int subject[observation_count];   // subject index
    int fragment[observation_count];  // fragment index
    int iteration_count;
    real subject_slope_delta_mean_estimate;
    real subject_slope_delta_variance_estimate;
    real baseline_slope;
    real fragment_slope_delta_reparameterized;
    real time_since_infection_variance_estimate;
} parameters {
} model {
}
generated quantities {
    vector[observation_count] time_since_infection;
    vector[subject_count] subject_slope_delta;
    vector[observation_count] total_slope;
    vector[observation_count] predicted_time_since_infection;
    for (subj in 1:subject_count){
        subject_slope_delta[subj] = normal_rng( subject_slope_delta_mean_estimate , subject_slope_delta_variance_estimate );
        for (i in 1:observation_count){
              total_slope[i] = (baseline_slope + subject_slope_delta[subj] + fragment_slope_delta_reparameterized);
              predicted_time_since_infection[i] = total_slope[i] * apd[i]; // linear function relating apd to time_since_infection
              time_since_infection[i] = normal_rng(predicted_time_since_infection[i], time_since_infection_variance_estimate);
        }
    }
}

