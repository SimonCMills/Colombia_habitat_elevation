// hierarchical habitat-detection model for elevational associations
// model: within-species asymmetric habitat scaling effect between upper and lower 
// limits, and habitat effect interacts with elevational range position

functions{
    real partial_sum(int[,] det_slice, 
                     int start, int end, 
                     int[] n_visit,
                     vector b0, 
                     vector b00,
                     vector b1,
                     vector b2,
                     vector b3,
                     vector b4,
                     vector d0, 
                     vector d1,
                     vector d2,
                     vector rMid, 
                     vector rScale,
                     row_vector[] vis_cov,
                     int[] id_sp, 
                     int[] id_site_sp,
                     int[,] id_obsvr,
                     int[] Q, 
                     vector ele,
                     vector habitat, 
                     vector midpoint) {
        // indexing vars
        int len = end - start + 1;
        int r0 = start - 1;
        
        vector[len] lp;
        real logit_psi;
        real eff_ele;
        for (r in 1:len) {
            int n_visit_r = n_visit[r0+r]; // number of visits to point
            row_vector[n_visit_r] logit_theta; 
            // elevational distribution component
            if(ele[r0+r] < rMid[id_sp[r0+r]])    
                eff_ele = b0[id_sp[r0+r]]  +
                    exp(rScale[id_sp[r0+r]] - b2[id_sp[r0+r]]*habitat[r0+r]) * 
                    -(ele[r0+r]-rMid[id_sp[r0+r]])^2 + 
                    b4[id_sp[r0+r]]*midpoint[r0+r];
            else 
                eff_ele = b0[id_sp[r0+r]]  +
                    exp(rScale[id_sp[r0+r]] - b3[id_sp[r0+r]]*habitat[r0+r]) * 
                    -(ele[r0+r]-rMid[id_sp[r0+r]])^2 + 
                    b4[id_sp[r0+r]]*midpoint[r0+r];
            // elevation + species:cluster ranef
            logit_psi = b00[id_site_sp[r0+r]] + eff_ele;
            // detection: species + time of day + observer
            logit_theta = d0[id_sp[r0+r]] + d1[id_sp[r0+r]]*vis_cov[r0+r, 1:n_visit_r] + 
                to_row_vector(d2[id_obsvr[r0+r, 1:n_visit_r]]);
            
            // likelihood
            if (Q[r0 + r] == 1) 
                lp[r] = log_inv_logit(logit_psi) +
                    bernoulli_logit_lpmf(det_slice[r, 1:n_visit_r] | logit_theta);
            else lp[r] = log_sum_exp(
                log_inv_logit(logit_psi) + sum(log1m_inv_logit(logit_theta)),
                log1m_inv_logit(logit_psi));
        } 
        return sum(lp);
    }
}

data {
    // dimensions
    int<lower=1> n_species; //number of species
    int<lower=1> n_sites_species; //number of species:site combinations
    int<lower=1> n_points; //number of species
    int<lower=1> n_tot; // nrows in df
    int<lower=1> n_visit[n_tot]; //variable number of visits
    // note: 4 is the maximum number of visits (visit-specific matrices all have 
    // this maximum dimension)
    
    // indexing variables
    int<lower=1> id_sp[n_tot]; // species ID
    int<lower=1> id_site_sp[n_tot]; // site:species ID
    int id_obsvr[n_tot, 4]; // observer
    int<lower=0, upper=1> Q[n_tot]; // detection/nondetection across visits
    
    // data & covariates
    row_vector[4] vis_cov[n_tot]; // time of day
    int det_data[n_tot, 4]; // detection history
    vector[n_tot] ele; //scaled elevations 
    vector[n_tot] habitat;
    vector[n_tot] midpoint; //range midpoint (a priori)
    int<lower=1> grainsize;
}
parameters {
    // occupancy: species intercept and site:species intercept
    real mu_b0;
    real<lower=0> sigma_b0;
    vector[n_species] b0_raw;
    
    real<lower=0> sigma_b00;
    vector[n_sites_species] b00_raw;
    
    // occupancy: species' habitat intercept and range-scaling effects
    real mu_b1;
    real<lower=0> sigma_b1;
    vector[n_species] b1_raw;

    real mu_b2;
    real<lower=0> sigma_b2;
    vector[n_species] b2_raw;
    
    real mu_b3;
    real<lower=0> sigma_b3;
    vector[n_species] b3_raw;
    
    real mu_b4;
    real<lower=0> sigma_b4;
    vector[n_species] b4_raw;
    
    // detection: species' intercept and time-of-day effect
    real mu_d0;
    real<lower=0> sigma_d0;
    vector[n_species] d0_raw;
    
    real mu_d1;
    real<lower=0> sigma_d1;
    vector[n_species] d1_raw;
    
    vector[3] d2; //observer effect
    
    // occupancy: range position (scaling and centre)
    real mu_rMid;
    real<lower=0> sigma_rMid;
    vector[n_species] rMid_raw;
    
    real mu_rScale;
    real<lower=0> sigma_rScale;
    vector[n_species] rScale_raw;
}
transformed parameters{
    // occupancy: intercept-terms
    vector[n_species] b0 = mu_b0 + b0_raw * sigma_b0;
    vector[n_sites_species] b00 = b00_raw * sigma_b00;
    // occupancy: habitat & range scaling
    vector[n_species] b1 = mu_b1 + b1_raw * sigma_b1;
    vector[n_species] b2 = mu_b2 + b2_raw * sigma_b2;
    vector[n_species] b3 = mu_b3 + b3_raw * sigma_b3;
    vector[n_species] b4 = mu_b4 + b4_raw * sigma_b4;
    vector[n_species] rMid = mu_rMid + rMid_raw * sigma_rMid;
    vector[n_species] rScale = mu_rScale + rScale_raw * sigma_rScale;
    // detection: varying species intercept
    vector[n_species] d0 = mu_d0 + d0_raw * sigma_d0;
    // detection: varying species time of day slope
    vector[n_species] d1 = mu_d1 + d1_raw * sigma_d1;
}
model {
    //Likelihood
    target += reduce_sum_static(partial_sum, det_data, grainsize, n_visit, b0, b00, 
        b1, b2, b3, b4, d0, d1, d2, rMid, rScale, vis_cov, id_sp, id_site_sp, id_obsvr, 
        Q, ele, habitat, midpoint);
    
    // Hyper-priors
    mu_b0 ~ normal(0,1.5);
    mu_b1 ~ normal(0,1.5);
    mu_b2 ~ normal(0,1.5);
    mu_b3 ~ normal(0,1.5);
    mu_b4 ~ normal(0,1.5);
    mu_d0 ~ normal(0,1.5);
    mu_d1 ~ normal(0,1.5);
    mu_rMid ~ normal(0,2);
    mu_rScale ~ normal(0,2);
    
    sigma_b0 ~ normal(0,1.5);
    sigma_b00 ~ normal(0,1.5);
    sigma_b1 ~ normal(0,1.5);
    sigma_b2 ~ normal(0,1.5);
    sigma_b3 ~ normal(0,1.5);
    sigma_b4 ~ normal(0,1.5);
    sigma_d0 ~ normal(0,1.5);
    sigma_d1 ~ normal(0,1.5);
    sigma_rMid ~ normal(0, 2);
    sigma_rScale ~ normal(0, 2);
    
    //Random Effects
    b0_raw ~ normal(0, 1);
    b00_raw ~ normal(0, 1);
    b1_raw ~ normal(0, 1);
    b2_raw ~ normal(0, 1);
    b3_raw ~ normal(0, 1);
    b4_raw ~ normal(0, 1);
    d0_raw ~ normal(0, 1);
    d1_raw ~ normal(0, 1);
    rMid_raw ~ normal(0, 1);
    rScale_raw ~ normal(0,1);
    
    // Observer effect
    d2 ~ normal(0, 1.5);
}
