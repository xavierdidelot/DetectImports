functions {
    //eigenfunction
    vector precompute_basis(vector T, real L, int j) {
        vector[size(T)] basis;
        real lambda = j*pi()/(2*L);
        basis = 1/sqrt(L) * sin(lambda*(T+L));
        return (basis);
    }

    //Spectral density function associated with Matern kernel with nu=3/2, dimension D=1 and Euclidiean distance
    real spec_dens_matern(real x, real alpha, real l) {
        real dens = 4 * (alpha^2) * (sqrt(3)/l)^3 * 1/((sqrt(3)/l)^2 + x^2)^2;
        return(dens);
    }

    real sqrt_spec_dens_matern(real x, real alpha, real l) {
        real dens = 2 * alpha * (sqrt(3)/l)^(1.5) * 1/((sqrt(3)/l)^2 + x^2);
        return(dens);
    }

    vector vec_sqrt_spd_matern(real rho, real L, int K) {
         return 2 * ((sqrt(3)/rho)^1.5) * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(K, 1, K))^2);
    }
}

data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    vector[N] T_s; // sampling times
    int<lower=1> M;// number of basis functions
    real<lower=1> c;// boundary
}

transformed data {
    // transform data to have domain[-1,1]
    vector[N] T_centered = (to_vector(T_s) - min(T_s));
    T_centered = T_centered - (max(T_centered)/2);
    real S=max(T_centered);
    T_centered = T_centered / S;
    real L = c;
    matrix[N, M] basis;
    for (idx in 1:M) {
        basis[:, idx] = precompute_basis(T_centered, L, idx);
    }
}

parameters {
    vector[M] f_tilde; // weights
    real<lower = 0> alpha;  // kernel scale (ie marginal std)
    real<lower = 0> l; // kernel length scale
}

transformed parameters {
    vector[N] a_coeffs;
    vector[N] coal_means;
    vector[M] spec_dens;

    spec_dens = vec_sqrt_spd_matern(l, L, M);
    a_coeffs = basis*(alpha*spec_dens.*f_tilde);
    coal_means = exp(a_coeffs);
}

model {
    alpha ~ normal(0,5);
    l ~ inv_gamma(5,5);
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_means);
}
