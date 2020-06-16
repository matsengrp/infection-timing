data{
    int<lower=1> observation_count; // number of observations
    int<lower=1> subject_count;   // number of subjects
    int<lower=1> fragment_count;  // number of fragments
    real observed_time[observation_count];   // measurement of outcome variable
    real apd[observation_count];    // predictor variable
    int subject[observation_count];   // subject index
    int fragment[observation_count];  // fragment index
}
parameters{
    real<lower=0> time_since_infection_variance_estimate;    // outcome time since infection variance (outcome uncertainty)
    real<lower=0> time_since_infection[observation_count];    // outcome
    real<lower=0, upper=150> baseline_slope;   // slope before subject/fragment slope changes
    vector[subject_count] subject_slope_delta;  // vector of subject specific slope changes
    real subject_slope_delta_mean_estimate;    // subject specific slope change mean
    real<lower=0> subject_slope_delta_variance_estimate; // subject specific slope change standard deviation
    vector[fragment_count] fragment_slope_delta;    // vector of fragment specific slope changes
    real fragment_slope_delta_mean_estimate;   // fragment specific slope change mean
    real<lower=0> fragment_slope_delta_variance_estimate;    // fragment specific slope change standard deviation
    vector<lower=0, upper=0.75>[subject_count] observed_time_to_time_since_infection_correction;    // conversion factor between observed time (age at sampling time) measurements and time since infection output
}
model{
    vector[observation_count] total_slope;
    vector[observation_count] predicted_time_since_infection;
    vector[observation_count] predicted_observed_time;
    fragment_slope_delta_mean_estimate ~ normal( 0 , 1 );
    fragment_slope_delta_variance_estimate ~ cauchy( 0 , 5 );
    subject_slope_delta_mean_estimate ~ normal( 0 , 1 );
    subject_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    fragment_slope_delta ~ normal( fragment_slope_delta_mean_estimate , fragment_slope_delta_variance_estimate );
    subject_slope_delta ~ normal( subject_slope_delta_mean_estimate , subject_slope_delta_variance_estimate );
    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta[fragment[i]]);
    }
    for ( i in 1:observation_count ) {
        predicted_time_since_infection[i] = total_slope[i] * apd[i]; // linear function relating apd to time_since_infection
    }
    observed_time_to_time_since_infection_correction ~ normal(0.5, 0.25);
    time_since_infection_variance_estimate ~ cauchy( 0 , 0.5 );
    time_since_infection ~ normal( predicted_time_since_infection , time_since_infection_variance_estimate);
    for ( i in 1:observation_count ) {
        predicted_observed_time[i] = time_since_infection[i] - observed_time_to_time_since_infection_correction[subject[i]]; // Conversion between predicted_observed_time (predicted age at sampling time) and time_since_infection. Here predicted_observed_time is calculated from the difference between time_since_infection and a subject specific observed_time_to_time_since_infection_correction
    }
    observed_time ~ normal(predicted_observed_time, 0.1);   // the observed_time (actual measured age at sampling time) is modeled as the predicted_observed_time value with some noise
}
generated quantities{
    vector[observation_count] total_slope;
    vector[observation_count] predicted_time_since_infection;
    vector[observation_count] predicted_observed_time;
    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta[fragment[i]]);
    }
    for ( i in 1:observation_count ) {
        predicted_time_since_infection[i] = total_slope[i] * apd[i];
    }
    for ( i in 1:observation_count ) {
        predicted_observed_time[i] = time_since_infection[i] - observed_time_to_time_since_infection_correction[subject[i]];
    }
}


