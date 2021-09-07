data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    real<lower = 0> T_s[N]; // sampling times
    real log_scale;
}

parameters {
    real<lower = 0> alpha;
    real log_mean;
    real f_tilde;
}
transformed parameters {
    real log_coal_mean = f_tilde*alpha+log_mean;
    real coal_mean = exp(log_coal_mean);
}

model {
    log_mean ~ normal(log_scale, log_scale);
    alpha ~ gamma(2,2);
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_mean);
}