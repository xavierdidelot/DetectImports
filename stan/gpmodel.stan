functions {
   vector precompute_basis(vector T, real b, int j) {
       vector[size(T)] basis;
       real lambda = j*pi()/(2*b);
       basis = 1/sqrt(b) * sin(lambda*(T+b));
       return (basis);
   }
   real spec_dens_rbf(real x, real alpha, real l) {
       return(alpha*sqrt(2*pi())*l*exp(-(l^2)*(x^2)/2));
   }
   real spec_dens_matern(real x, real alpha, real l) {
       real dens = alpha*(2*sqrt(pi())*3^(3/2))/((1/2)*sqrt(pi())*l^3)*((3/(l^2))+x^2)^(-2);
       return(dens);
   }
}

data {
    int<lower = 1> N; // sample number
    real<lower = 0> intervals[N]; // coalescent intervals
    real<lower = 0> T_s[N]; // sampling times
    real<lower = 0> shape;
    real<lower = 0> scale;
    int<lower=1> M;
    real<lower=1> c;
}

transformed data {
    vector[N] T_centered = (to_vector(T_s) - min(T_s));
    #T_centered = T_centered - (max(T_centered)/2);
    T_centered = 2*T_centered/(max(T_centered)-min(T_centered))-1;
    real L = c*max(T_centered);
    matrix[N, M] basis;
    for (idx in 1:M) {
        basis[:, idx] = precompute_basis(T_centered, L ,idx);
    }
}

parameters {
    vector[M] f_tilde;
    real<lower = 0> alpha;  // kernel sigma
    real<lower = 0.3> l; // kernel length scale
}

transformed parameters {
    vector[N] a_coeffs;
    vector[N] coal_means;
    vector[M] spec_dens;
    {
        for(idx in 1:M) {
            spec_dens[idx] = sqrt(spec_dens_matern(idx*pi()/(2*L), alpha, l));
        }
        //print(spec_dens);
        a_coeffs = (basis)*(spec_dens.*f_tilde); 
    }
    coal_means = exp(a_coeffs);
}

model {
    alpha ~ cauchy(0,1);
    l ~ inv_gamma(shape,scale);
    f_tilde ~ normal(0,1);
    intervals ~ exponential(coal_means);
}

generated quantities {
    vector[N] e_tilde;
    real f[N];
    {
        e_tilde = to_vector(intervals)./coal_means;
        f = exponential_rng(coal_means);
    }
}

