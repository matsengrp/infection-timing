data{
    int<lower=1> observation_count; // number of observations
    int<lower=1> subject_count;   // number of subjects
    int<lower=1> fragment_count;  // number of fragments
    real observed_time_since_infection[observation_count];   // measurement of outcome variable
    real apd[observation_count];    // predictor variable
    int subject[observation_count];   // subject index
    int fragment[observation_count];  // fragment index
    int is_utero[observation_count];
    int is_post[observation_count];
    real vload[observation_count];
}
parameters{
    real<lower=0> time_since_infection_variance_estimate;    // outcome time since infection variance (outcome uncertainty)
    real<lower=0> time_since_infection[observation_count];    // outcome
    real<lower=0, upper=150> baseline_slope;   // slope before subject/fragment slope changes
    real vload_slope_delta;   // slope before subject/fragment slope changes
    real vload_slope_delta_mean_estimate;
    real<lower=0> vload_slope_delta_variance_estimate;
    vector[subject_count] subject_slope_delta;  // vector of subject specific slope changes
    real subject_slope_delta_mean_estimate;    // subject specific slope change mean
    real<lower=0> subject_slope_delta_variance_estimate; // subject specific slope change standard deviation
    vector[fragment_count] fragment_slope_delta_raw;    // vector of fragment specific slope changes
    real fragment_slope_delta_mean_estimate;   // fragment specific slope change mean
    real<lower=0> fragment_slope_delta_variance_estimate;    // fragment specific slope change standard deviation
    vector<lower=0, upper=0.75>[subject_count] observed_time_to_time_since_infection_correction_utero;    // conversion factor between observed time (age at sampling time) measurements and time since infection output
}
transformed parameters{
    vector[fragment_count] fragment_slope_delta_reparameterized;
    fragment_slope_delta_reparameterized = fragment_slope_delta_mean_estimate + fragment_slope_delta_variance_estimate * fragment_slope_delta_raw;
}
model{
    vector[observation_count] total_slope;
    vector[observation_count] predicted_time_since_infection;
    vector[observation_count] predicted_observed_time_since_infection;
    vector[observation_count] observed_time_to_time_since_infection_correction;
    fragment_slope_delta_mean_estimate ~ normal( 0 , 1 );
    fragment_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    fragment_slope_delta_raw ~ normal(0, 1);
    subject_slope_delta_mean_estimate ~ normal( 0 , 1 );
    subject_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    subject_slope_delta ~ normal( subject_slope_delta_mean_estimate , subject_slope_delta_variance_estimate );
    vload_slope_delta_mean_estimate ~ normal( 0 , 1 );
    vload_slope_delta_variance_estimate ~ cauchy( 0 , 20 );
    vload_slope_delta ~ normal( vload_slope_delta_mean_estimate , vload_slope_delta_variance_estimate );

    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta_reparameterized[fragment[i]]);
    }
    for ( i in 1:observation_count ) {
        predicted_time_since_infection[i] = total_slope[i] * apd[i] + vload_slope_delta * apd[i]*vload[i]; // linear function relating apd to time_since_infection
    }
    observed_time_to_time_since_infection_correction_utero ~ beta(1,5); // subjects infected in utero were most likely infected in the third trimester (of pregnancy)

    time_since_infection_variance_estimate ~ cauchy( 0 , 0.5 );
    time_since_infection ~ normal( predicted_time_since_infection , time_since_infection_variance_estimate);
    for ( i in 1:observation_count ) {
        if ( is_utero[i] ){
          observed_time_to_time_since_infection_correction[subject[i]] = observed_time_to_time_since_infection_correction_utero[subject[i]];
        } else if ( is_post[i] ){
          observed_time_to_time_since_infection_correction[subject[i]] = 0;
        }
    }
    for ( i in 1:observation_count ) {
        predicted_observed_time_since_infection[i] = time_since_infection[i] - observed_time_to_time_since_infection_correction[subject[i]];
    }
    observed_time_since_infection ~ normal(predicted_observed_time_since_infection, 0.1);   // the observed_time (actual measured age at sampling time) is modeled as the predicted_observed_time value with some noise
}
generated quantities{
    vector[observation_count] total_slope;
    vector[observation_count] predicted_time_since_infection;
    vector[subject_count] observed_time_to_time_since_infection_correction;
    for ( i in 1:observation_count ) {
        total_slope[i] = (baseline_slope + subject_slope_delta[subject[i]] + fragment_slope_delta_reparameterized[fragment[i]]);
    }
    for ( i in 1:observation_count ) {
        predicted_time_since_infection[i] = total_slope[i] * apd[i] + vload_slope_delta * apd[i] * vload[i];
    }
    for ( i in 1:observation_count ) {
        if ( is_utero[i] ){
          observed_time_to_time_since_infection_correction[subject[i]] = observed_time_to_time_since_infection_correction_utero[subject[i]];
        } else if ( is_post[i] ){
          observed_time_to_time_since_infection_correction[subject[i]] = 0;
        }
    }
}
