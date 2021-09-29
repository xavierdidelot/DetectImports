functions {
   vector precompute_basis(vector T, real L, int j) {
       vector[size(T)] basis;
       real lambda = j*pi()/(2*L);
       basis = 1/sqrt(L) * sin(lambda*(T+L));
       return (basis);
   }
   real spec_dens_matern(real x, real alpha, real l) {
       real dens = 4 * (alpha^2) * (sqrt(3)/l)^3 * 1/((sqrt(3)/l)^2 + x^2)^2;
       return(dens);
   }
}

data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    vector[N] T_s; // sampling times
    real<lower = 0> shape;
    real<lower = 0> scale;
    int<lower=1> M;
    real<lower=1> c;
}

transformed data {
    vector[N] T_centered = (to_vector(T_s) - min(T_s));
    T_centered = T_centered - (max(T_centered)/2);
    real S=max(T_centered);
    real L = c*S;
    matrix[N, M] basis;
    for (idx in 1:M) {
        basis[:, idx] = precompute_basis(T_centered, L ,idx);
    }
}

parameters {
    vector[M] f_tilde;
    real<lower = 0> alpha;  // kernel sigma
    real<lower = 0.01> l; // kernel length scale
}

transformed parameters {
    vector[N] a_coeffs;
    vector[N] coal_means;
    vector[M] spec_dens;
    {
        for(idx in 1:M) {
            spec_dens[idx] = sqrt(spec_dens_matern(idx*pi()/(2*L), alpha, l));
        }
        a_coeffs = (basis)*(spec_dens.*f_tilde);
    }
    coal_means = exp(a_coeffs);
}

model {
    alpha ~ gamma(2,2);
    l ~ inv_gamma(shape,scale);
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_means);
}
