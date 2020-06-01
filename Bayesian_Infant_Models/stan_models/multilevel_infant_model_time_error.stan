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
    vector[N_subject_index] subject_slope;  // vector of subject specific slope changes
    vector[N_fragment_index] fragment_slope;    // vector of fragment specific slope changes
    vector<lower=0,upper=0.75>[N_subject_index] time_correction;    // conversion factor between time(age) measurements and time since infection output
    real<lower=0> fragment_slope_sd;    // fragment specific slope standard deviation
    real fragment_slope_mean;   // fragment specific slope mean
    real<lower=0> subject_slope_sd; // subject specific slope standard deviation
    real subject_slope_mean;    // subject specific slope mean
    real<lower=0,upper=150> population_avg_slope;   // average slope
}
model{
    vector[N] total_slope;
    vector[N] population_mean;
    vector[N] time_adjustment;
    fragment_slope_sd ~ cauchy( 0 , 5 );
    fragment_slope_mean ~ normal( 0 , 1 );
    subject_slope_sd ~ cauchy( 0 , 20 );
    subject_slope_mean ~ normal( 0 , 1 );
    fragment_slope ~ normal( fragment_slope_mean , fragment_slope_sd );
    subject_slope ~ normal( subject_slope_mean , subject_slope_sd );
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[fragment_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_slope[i] * apd[i];
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time_since_infection ~ normal( population_mean , population_sd );
    for ( i in 1:N ) {
        time_adjustment[i] = time_since_infection[i] - time_correction[subject_index[i]]; // Conversion between age at sampling time (time_adjustment) and time since infection. Here time_adjustment is the predicted time (age) calculated from the time_since_infection minus a subject specific time_correction
    }
    time ~ normal(time_adjustment, 0.01);   // the measured age at sampling time is modeled as the time_adjustment value (time (age) estimate) with some noise...
}
generated quantities{
    vector[N] total_slope;
    vector[N] population_mean;
    vector[N] time_adjustment;
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[fragment_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_slope[i] * apd[i];
    }
    for ( i in 1:N ) {
        time_adjustment[i] = time_since_infection[i] - time_correction[subject_index[i]]; //This means that we aren't directly modelling time_since_infection
    }
}


//This model does not have any divergent iterations, however, it claims to have one Rhat value of NA (at a data point where APD  0)


