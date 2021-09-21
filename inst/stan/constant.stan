data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    real<lower = 0> T_s[N]; // sampling times
}

parameters {
    real<lower = 0> coal_mean;
}

model {
    coal_mean ~ inv_gamma(0.001,0.001);
    intervals ~ exponential(coal_mean);
}
