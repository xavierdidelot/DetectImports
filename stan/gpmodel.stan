data {
    int<lower = 1> N; // sample number
    real<lower = 0> sigma;  // kernel sigma
    real<lower = 0> l; // kernel length scale
    real<lower = 0> intervals[N]; // coalescent intervals
    real<lower = 0> T_s[N]; // sampling times
}

transformed data {
    matrix[N, N] cov = cov_exp_quad(T_s, sigma, l) + diag_matrix(rep_vector(1e-10, N));
    matrix[N, N] L_cov = cholesky_decompose(cov);
    
}

parameters {
    vector[N] f_tilde;
}

transformed parameters {
    vector[N] log_mean_c = L_cov*f_tilde; 
    vector[N] coal_mean = exp(log_mean_c);
}


model {
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_mean);
}

generated quantities {
    real f[N] = exponential_rng(coal_mean);
}