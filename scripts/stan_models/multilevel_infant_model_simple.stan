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
    real<lower=0, upper=150> baseline_slope;   // slope before subject/fragment slope changes
    real<lower=0> time_variance;
}
model{
    for ( i in 1:observation_count ) {
        observed_time[i] ~ normal(baseline_slope * apd[i], time_variance); 
    }
}


