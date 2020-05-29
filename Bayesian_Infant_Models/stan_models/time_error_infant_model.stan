data{
    int<lower=1> N; // number of observations
    int<lower=1> N_subject_index;   // number of subject clusters
    int<lower=1> N_frag_index;  // number of fragment clusters
    real time[N];   // measurement of time 
    real apd[N];    // measurement of apd
    int subject_index[N];   
    int frag_index[N];
}
parameters{
    vector[150] time_since_infection;
    real<lower=0> population_sd;
    real<lower=0,upper=0.75> time_error;
    vector[N_subject_index] subject_slope;
    vector[N_frag_index] fragment_slope;
    real subject_slope_mean;
    real<lower=0> subject_slope_sd;
    real<lower=0,upper=150> population_avg_slope;
    real fragment_slope_mean;
    real<lower=0> fragment_slope_sd;
}
model{
    vector[N] total_slope;
    vector[N] population_mean;
    fragment_slope_sd ~ cauchy( 0 , 5 );
    fragment_slope_mean ~ normal( 0 , 1 );
    // population_avg_slope ~ uniform( 0 , 150 );
    subject_slope_sd ~ cauchy( 0 , 20 );
    subject_slope_mean ~ normal( 0 , 1 );
    fragment_slope ~ normal( fragment_slope_mean , fragment_slope_sd );
    subject_slope ~ normal( subject_slope_mean , subject_slope_sd );
    // time_error ~ uniform( 0 , 0.75 );
    time ~ cauchy( time_since_infection - time_error , 0.05 );
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_slope[i] * apd[i];
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time_since_infection ~ normal( population_mean , population_sd );
}
generated quantities{
    vector[N] total_slope;
    vector[N] population_mean;
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_slope[i] * apd[i];
    }
}