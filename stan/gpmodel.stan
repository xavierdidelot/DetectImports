
data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    real<lower = 0> T_s[N]; // sampling times
    real<lower = 0> shape;
    real<lower = 0> scale;
    int<lower=1 M> order;
}

functions {
   vector precompute_basis(real[] Ts, real L, int j) {
       vector[size(Ts)] basis;
       real lambda = (j*pi/(2*L))^2;
       basis = 1/sqrt(L) * sin(sqrt(lambda)*(Ts+L));
       return (basis);
   }
   real spec_dens_rbf(real x, real alpha, real l) {
       return(alpha*sqrt(2*pi)*l*exp(-(l^2)*(x^2)/2));
   }
   real spec_dens_matern(real x, real alpha, real l) {
       real dens = alpha*(2*sqrt(pi)*tgamma(2)*3^(3/2))/(1/2*sqrt(pi)*l^3)*(3/l^2+x^2)^(-2);
       return(dens);
   }
}

parameters {
    vector[N] f_tilde;
    real<lower = 0> sigma;  // kernel sigma
    real<lower = 0> l; // kernel length scale
}

transformed parameters {
    vector[N] a_coeffs;
    vector[N] coal_means;
    {
        matrix[N, N] cov = cov_exp_quad(T_s, sigma, l) + diag_matrix(rep_vector(1e-8, N));
        matrix[N, N] L_cov = cholesky_decompose(cov);
        a_coeffs = L_cov*f_tilde; 
    }
    coal_means = exp(a_coeffs);
}

model {
    sigma ~ cauchy(0,1);
    l ~ inv_gamma(shape,scale);
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_means);
}

generated quantities {
    vector[N] e_tilde;
    {
        e_tilde = to_vector(intervals)./coal_means;
    }
}

