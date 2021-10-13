data{
    int<lower=1> observation_count; // number of observations
    int<lower=1> subject_count;   // number of subjects
    int<lower=1> fragment_count;  // number of fragments
    real observed_time_since_infection[observation_count];   // measurement of outcome variable
    real apd[observation_count];    // predictor variable
    int subject[observation_count];   // subject index
    int fragment[observation_count];  // fragment index
}
parameters{
    real<lower=0, upper=150> baseline_slope;   // slope before subject/fragment slope changes
    vector[subject_count] subject_slope_delta;  // vector of subject specific slope changes
    real subject_slope_delta_mean_estimate;    // subject specific slope change mean
    real<lower=0> subject_slope_delta_variance_estimate; // subject specific slope change standard deviation
    vector[fragment_count] fragment_slope_delta_raw;    // vector of fragment specific slope changes
    real fragment_slope_delta_mean_estimate;   // fragment specific slope change mean
    real<lower=0> fragment_slope_delta_variance_estimate;    // fragment specific slope change standard deviation
    real<lower=0> time_variance;
}
transformed parameters{
    vector[fragment_count] fragment_slope_delta_reparameterized;
    fragment_slope_delta_reparameterized = fragment_slope_delta_mean_estimate + fragment_slope_delta_variance_estimate * fragment_slope_delta_raw;
}
model{
    vector[observation_count] total_slope;
    fragment_slope_delta_mean_estimate ~ normal( 0 , 1 );
    fragment_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    fragment_slope_delta_raw ~ normal(0, 1);
    subject_slope_delta_mean_estimate ~ normal( 0 , 1 );
    subject_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    subject_slope_delta ~ normal( subject_slope_delta_mean_estimate , subject_slope_delta_variance_estimate );
    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta_reparameterized[fragment[i]]);
    }
    for ( i in 1:observation_count ) {
        observed_time_since_infection[i] ~ normal(total_slope[i] * apd[i], time_variance); 
    }
}
generated quantities{
    vector[observation_count] total_slope;
    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta_reparameterized[fragment[i]]);
    }
}

