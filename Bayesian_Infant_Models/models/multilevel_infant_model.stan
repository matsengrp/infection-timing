data{
    // Define variables in data
    // Number of observations
    int<lower=1> N;
    // Number of Subject IDs
    int<lower=1> N_subject_index;
    // Number of Fragment IDs
    int<lower=1> N_frag_index;
    // Continuous outcome
    real time[N];
    // Continuous predictor
    real apd[N];
    // Subject IDs
    int subject_index[N];
    // Fragment IDs
    int frag_index[N];
}
parameters{
    // Define parameters to estimate
    // Population standard deviation
    real<lower=0> population_sd;
    // Population total intercept
    real<lower=-0.75,upper=0> total_intercept;
    // Varying slope subject
    vector[N_subject_index] subject_slope;
    // Varying slope fragment
    vector[N_frag_index] fragment_slope;
    // Subject specific slope mean
    real subject_slope_mean;
    // Subject specific slope sd
    real<lower=0> subject_slope_sd;
    // Population slope mean
    real<lower=0,upper=150> population_avg_slope;
    // Fragment specific slope mean
    real fragment_slope_mean;
    // Fragment specific slope sd
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
    // total_intercept ~ uniform( -0.75 , 0 );
    for ( i in 1:N ) {
        total_slope[i] = population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]];
    }
    for ( i in 1:N ) {
        population_mean[i] = total_intercept + total_slope[i] * apd[i];
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time ~ normal( population_mean , population_sd );
}
generated quantities{
    vector[N] total_slope;
    vector[N] population_mean;
    for ( i in 1:N ) {
        total_slope[i] = population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]];
    }
    for ( i in 1:N ) {
        population_mean[i] = total_intercept + total_slope[i] * apd[i];
    }
}


data{
    // Define variables in data
    // Number of observations
    int<lower=1> N;
    // Number of Subject IDs
    int<lower=1> N_subject_index;
    // Number of Fragment IDs
    int<lower=1> N_frag_index;
    // Continuous outcome
    real time[N];
    // Continuous predictor
    real apd[N];
    // Subject IDs
    int subject_index[N];
    // Fragment IDs
    int frag_index[N];
}
parameters{
    // Define parameters to estimate
    // Population standard deviation
    real<lower=0> population_sd;
    // Varying slope subject
    vector[N_subject_index] subject_slope;
    // Varying intercept subject
    vector<lower=-0.75,upper=0>[N_subject_index] subject_intercept;
    // Varying slope fragment
    vector[N_frag_index] fragment_slope;
    // Subject specific slope mean
    real subject_slope_mean;
    // Subject specific intercept mean
    real subject_intercept_mean;
    // Subject correlation matrix
    corr_matrix[2] subject_correlation_matrix;
    // Subject specific sd
    vector<lower=0>[2] subject_sd;
    // Population slope mean
    real<lower=0,upper=150> population_avg_slope;
    // Fragment specific slope mean
    real fragment_slope_mean;
    // Fragment specific slope sd
    real<lower=0> fragment_slope_sd;
}
transformed parameters{
    vector[2] v_subject_interceptsubject_slope[N_subject_index];
    vector[2] Mu_subject_intercept_meansubject_slope_mean;
    cov_matrix[2] SRS_subject_sdsubject_correlation_matrix;
    for ( j in 1:N_subject_index ) {
        v_subject_interceptsubject_slope[j,1] = subject_intercept[j];
        v_subject_interceptsubject_slope[j,2] = subject_slope[j];
    }
    for ( j in 1:2 ) {
        Mu_subject_intercept_meansubject_slope_mean[1] = subject_intercept_mean;
        Mu_subject_intercept_meansubject_slope_mean[2] = subject_slope_mean;
    }
    SRS_subject_sdsubject_correlation_matrix = quad_form_diag(subject_correlation_matrix,subject_sd);
}
model{
    vector[N] total_slope;
    vector[N] population_mean;
    fragment_slope_sd ~ cauchy( 0 , 5 );
    fragment_slope_mean ~ normal( 0 , 1 );
    // population_avg_slope ~ uniform( 0 , 150 )T[0, ];
    subject_sd ~ cauchy( 0 , 20 );
    subject_correlation_matrix ~ lkj_corr( 1 );
    subject_intercept_mean ~ normal( 0 , 1 )T[-0.75, 0];
    subject_slope_mean ~ normal( 0 , 1 );
    fragment_slope ~ normal( fragment_slope_mean , fragment_slope_sd );
    v_subject_interceptsubject_slope ~ multi_normal( Mu_subject_intercept_meansubject_slope_mean , SRS_subject_sdsubject_correlation_matrix );
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = subject_intercept[subject_index[i]] + total_slope[i] * apd[i];
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time ~ normal( population_mean , population_sd );
}
generated quantities{
    vector<lower=0>[N] total_slope;
    vector[N] population_mean;
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = subject_intercept[subject_index[i]] + total_slope[i] * apd[i];
    }
}




data{
    int<lower=1> N;
    int<lower=1> N_frag_index;
    int<lower=1> N_subject_index;
    real time[N];
    real apd[N];
    int subject_index[N];
    int frag_index[N];
}
parameters{
    real<lower=0> population_sd;
    vector[N_subject_index] subject_slope;
    vector[N_frag_index] fragment_slope;
    vector[N_subject_index] subject_intercept;
    real subject_slope_mean;
    real<lower=0> subject_slope_sd;
    real subject_intercept_mean;
    real<lower=0> subject_intercept_sd;
    real<lower=0,upper=150> population_avg_slope;
    real population_avg_int;
    real fragment_slope_mean;
    real<lower=0> fragment_slope_sd;
}
model{
    vector[N] total_intercept;
    vector[N] total_slope;
    vector[N] population_mean;
    fragment_slope_sd ~ cauchy( 0 , 5 );
    fragment_slope_mean ~ normal( 0 , 1 );
    population_avg_int ~ normal( 0 , 0.5 );
    // population_avg_slope ~ uniform( 0 , 150 );
    subject_intercept_sd ~ cauchy( 0 , 0.25 );
    subject_intercept_mean ~ normal( 0 , 1 );
    subject_slope_sd ~ cauchy( 0 , 20 );
    subject_slope_mean ~ normal( 0 , 1 );
    subject_intercept ~ normal( subject_intercept_mean , subject_intercept_sd );
    fragment_slope ~ normal( fragment_slope_mean , fragment_slope_sd );
    subject_slope ~ normal( subject_slope_mean , subject_slope_sd );
    for ( i in 1:N ) {
        total_intercept[i] = (population_avg_int + subject_intercept[subject_index[i]]);
    }
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_intercept[i] + total_slope[i] * apd[i];
    }
    population_sd ~ cauchy( 0 , 0.5 );
    time ~ normal( population_mean , population_sd );
}
generated quantities{
    vector[N] total_intercept;
    vector[N] total_slope;
    vector[N] population_mean;
    for ( i in 1:N ) {
        total_intercept[i] = (population_avg_int + subject_intercept[subject_index[i]]);
    }
    for ( i in 1:N ) {
        total_slope[i] = (population_avg_slope + subject_slope[subject_index[i]] + fragment_slope[frag_index[i]]);
    }
    for ( i in 1:N ) {
        population_mean[i] = total_intercept[i] + total_slope[i] * apd[i];
    }
}