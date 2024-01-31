data {
    int<lower = 1> N; // sample number
    array[N] real<lower = 0> intervals; // coalescent intervals
    array[N] real<lower = 0> T_s; // sampling times
}

parameters {
    real<lower = 0> coal_mean;
}

model {
    coal_mean ~ inv_gamma(0.001,0.001);
    intervals ~ exponential(coal_mean);
}
