// compare non-centered parameterisations for speed; compare speed for modelling
// covariance
// old-style specification, but define pars transform in the model block 

data {
    int<lower=1> n_pt; 
    vector[n_pt] x;
    int y[n_pt];
    int<lower=1> n_grp;
    int id_grp[n_pt];
}
parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_grp] b0_raw;
    
    real mu_b1;
    real<lower=0> sigma_b1;
    vector[n_grp] b1_raw;
}
model {
    vector[n_grp] b0 = mu_b0 + sigma_b0*b0_raw;
    vector[n_grp] b1 = mu_b1 + sigma_b1*b1_raw;
    for(i in 1:n_pt) {
        y[i] ~ bernoulli_logit(b0[id_grp[i]] + b1[id_grp[i]]*x[i]);
    }
    // priors
    b0_raw ~ std_normal();
    b1_raw ~ std_normal();
}
