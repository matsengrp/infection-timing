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
    vector[N_subject_index] subject_slope;  // vector of subject specific slope changes
    vector[N_fragment_index] fragment_slope;    // vector of fragment specific slope changes
    real<lower=0> fragment_slope_sd;    // fragment specific slope standard deviation
    real fragment_slope_mean;   // fragment specific slope mean
    real<lower=0> subject_slope_sd; // subject specific slope standard deviation
    real subject_slope_mean;    // subject specific slope mean
    real<lower=0,upper=150> population_avg_slope;   // average slope
    real error_intercept;   // average slope
}
model{
    vector[N] total_slope;
    vector[N] population_mean;
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
        population_mean[i] = total_slope[i] * apd[i] + error_intercept;
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time ~ normal( population_mean , population_sd );
}
generated quantities{
    vector[N] total_slope;
    vector[N] population_mean;
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[fragment_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_slope[i] * apd[i] + error_intercept;
    }
}



