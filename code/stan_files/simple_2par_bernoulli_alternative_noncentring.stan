// compare non-centered parameterisations for speed; compare speed for modelling
// covariance
// alternative non-centring

data {
    int<lower = 1> n_pt; 
    vector[n_pt] x;
    int y[n_pt];
    int<lower=1> n_grp;
    int id_grp[n_pt];
}
parameters {
    real mu_b0;
    real<lower=0> sigma_b0;
    vector<offset = mu_b0, multiplier = sigma_b0>[n_grp] b0;
    
    real mu_b1;
    real<lower=0> sigma_b1;
    vector<offset = mu_b1, multiplier = sigma_b1>[n_grp] b1;
}
model {
    for(i in 1:n_pt) {
        y[i] ~ bernoulli_logit(b0[id_grp[i]] + b1[id_grp[i]]*x[i]);
    }
    
    b0 ~ normal(mu_b0, sigma_b0);
    b1 ~ normal(mu_b1, sigma_b1);
}
