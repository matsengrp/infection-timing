data{
    int<lower=1> N; // number of observations
    vector[N] apd;    // predictor variable
    real time[N];
}
parameters{
    real time_since_infection[N];
    real<lower=0> total_slope;  // slope
    real<lower=0> time_noise;   // outcome time noise
    real<lower=0> time_correction_noise;
}
model {
    time ~ normal(time_since_infection, time_correction_noise);
    time_since_infection ~ normal(total_slope * apd, time_noise); 
    time_correction_noise ~ uniform(0,0.75);
    total_slope ~ uniform(0,300);
    time_noise ~ cauchy(0, 0.5);
}


