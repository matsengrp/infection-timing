data{
    int<lower=1> N; // number of observations
    int<lower=1> N_subject_index;   // number of subjects
    int<lower=1> N_fragment_index;  // number of fragments
    real time[N];   // predictor variable
    real apd[N];    // measurement of outcome variable
    int subject_index[N];   // subject index
    int fragment_index[N];  // fragment index
}
parameters{
    real<lower=0> population_sd;    // outcome time since infection noise (outcome uncertainty)
    real<lower=0> time_since_infection[N];    // outcome
    vector[N_subject_index] subject_shape_param;  // vector of subject specific slope changes
    vector[N_fragment_index] fragment_shape_param;    // vector of fragment specific slope changes
    vector<lower=0,upper=0.75>[N_subject_index] time_correction;    // conversion factor between time(age) measurements and time since infection output
    real<lower=0> fragment_shape_sd;    // fragment specific slope standard deviation
    real fragment_shape_mean;   // fragment specific slope mean
    real<lower=0> subject_shape_sd; // subject specific slope standard deviation
    real subject_shape_mean;    // subject specific slope mean
    real<lower=0> average_shape_parameter; // average shape parameter
    real<lower=0> up_shift;   // average shift upwards
    real<lower=0> right_shift;   // average shift right
}
transformed parameters {
    real shift_power_apd[N];

    for(i in 1:N) {
        shift_power_apd[i] = (apd[i] - right_shift)^3;
    }
}
model{
    vector[N] shape_param;
    vector[N] population_mean;
    vector[N] sampling_time_estimate;
    fragment_shape_sd ~ cauchy( 0 , 5 );
    fragment_shape_mean ~ normal( 0 , 100 );
    subject_shape_sd ~ cauchy( 0 , 5 );
    subject_shape_mean ~ normal( 0 , 100 );
    fragment_shape_param ~ normal( fragment_shape_mean , fragment_shape_sd );
    subject_shape_param ~ normal( subject_shape_mean , subject_shape_sd );
    for ( i in 1:N ) {
        shape_param[i] = (average_shape_parameter + subject_shape_param[subject_index[i]] + fragment_shape_param[fragment_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = shape_param[i]*shift_power_apd[i] + up_shift;
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time_since_infection ~ normal( population_mean , population_sd );
    for ( i in 1:N ) {
        sampling_time_estimate[i] = time_since_infection[i] - time_correction[subject_index[i]]; // Conversion between age at sampling time (sampling_time_estimate) and time since infection. Here sampling_time_estimate is the predicted time (age) calculated from the time_since_infection minus a subject specific time_correction
    }
    time ~ normal(sampling_time_estimate, 0.01);   // the measured age at sampling time is modeled as the sampling_time_estimate value (time (age) estimate) with some noise...
}
generated quantities{
    vector[N] shape_param;
    vector[N] population_mean;
    vector[N] sampling_time_estimate;
    for ( i in 1:N ) {
        shape_param[i] = (average_shape_parameter + subject_shape_param[subject_index[i]] + fragment_shape_param[fragment_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = shape_param[i]*shift_power_apd[i] + up_shift;
    }
    for ( i in 1:N ) {
        sampling_time_estimate[i] = time_since_infection[i] - time_correction[subject_index[i]]; //This means that we aren't directly modelling time_since_infection
    }
}
