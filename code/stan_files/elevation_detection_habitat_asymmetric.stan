// hierarchical habitat-detection model for elevational associations
// model 2: asymmetric habitat scaling effect between upper and lower limits
functions{
    real partial_sum(int[,] det_slice, 
                     int start, int end, 
                     int n_visit,
                     vector b0, 
                     vector b00,
                     vector b1,
                     vector b11,
                     vector b12,
                     vector d0, 
                     vector d1,
                     vector rMid, 
                     vector rScale,
                     row_vector[] vis_cov1, // time (10-stop block)
                     int[] id_sp, 
                     int[] id_site_sp,
                     int[] Q, 
                     vector ele,
                     vector habitat) {
        // indexing vars
        int len = end - start + 1;
        int r0 = start - 1;
        
        vector[len] lp;
        real logit_psi;
        row_vector[n_visit] logit_theta;
        real ele_scaled;
        for (r in 1:len) {
            // elevational distribution component (i.e. x-axis stuff)
            if(ele[r0+r] < rMid[id_sp[r0+r]])
                ele_scaled = exp(rScale[id_sp[r0+r]] - b11[id_sp[r0+r]]*habitat[r0+r]) * 
                    -(ele[r0+r]-rMid[id_sp[r0+r]])^2;
            else ele_scaled = exp(rScale[id_sp[r0+r]] - b12[id_sp[r0+r]]*habitat[r0+r]) * 
                -(ele[r0+r]-rMid[id_sp[r0+r]])^2;
            
            // add in y-axis stuff
            logit_psi = b0[id_sp[r0+r]] + b00[id_site_sp[r0+r]] + b1[id_sp[r0+r]]*habitat[r0+r] + ele_scaled;
            logit_theta = d0[id_sp[r0+r]] + d1[id_sp[r0+r]]*vis_cov1[r0+r];
            // likelihood
            if (Q[r0 + r] == 1) 
                lp[r] = log_inv_logit(logit_psi) +
                    bernoulli_logit_lpmf(det_slice[r] | logit_theta);
            else lp[r] = log_sum_exp(
                log_inv_logit(logit_psi) +
                    log1m_inv_logit(logit_theta[1]) +
                    log1m_inv_logit(logit_theta[2]) +
                    log1m_inv_logit(logit_theta[3]) +
                    log1m_inv_logit(logit_theta[4]),
                log1m_inv_logit(logit_psi));
        } 
        return sum(lp);
    }
}
data {
    // dimensions
    int<lower=1> n_visit; //fixed number of visits
    int<lower=1> n_species; //number of species
    int<lower=1> n_sites_species; //number of species:site combinations
    int<lower=1> n_points; //number of species
    int<lower=1> n_tot; // nrows in df
    
    // indexing variables
    int<lower=1> id_sp[n_tot]; // species ID
    int<lower=1> id_site_sp[n_tot]; // site:species ID
    int<lower=0, upper=1> Q[n_tot]; // detection/nondetection across visits
    
    // data & covariates
    row_vector[n_visit] vis_cov1[n_tot]; // time of day
    int<lower=0, upper=1> det_data[n_tot, n_visit]; // detection history
    vector[n_tot] ele; //scaled elevations 
    vector[n_tot] habitat;
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

    real mu_b11;
    real<lower=0> sigma_b11;
    vector[n_species] b11_raw;
    
    real mu_b12;
    real<lower=0> sigma_b12;
    vector[n_species] b12_raw;
    
    // detection: species' intercept and time-of-day effect
    real mu_d0;
    real<lower=0> sigma_d0;
    vector[n_species] d0_raw;
    
    real mu_d1;
    real<lower=0> sigma_d1;
    vector[n_species] d1_raw;
    
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
    vector[n_species] b11 = mu_b11 + b11_raw * sigma_b11;
    vector[n_species] b12 = mu_b12 + b12_raw * sigma_b12;
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
        b1, b11, b12, d0, d1, rMid, rScale, vis_cov1, id_sp, id_site_sp, Q, ele, habitat);
    
    // Hyper-priors:
    mu_b0 ~ normal(0,5);
    mu_b1 ~ normal(0,5);
    mu_b11 ~ normal(0,5);
    mu_b12 ~ normal(0,5);
    mu_d0 ~ normal(0,5);
    mu_d1 ~ normal(0,5);
    mu_rMid ~ normal(0,5);
    mu_rScale ~ normal(0,5);
    
    sigma_b0 ~ normal(0,5);
    sigma_b00 ~ normal(0,5);
    sigma_b1 ~ normal(0,5);
    sigma_b11 ~ normal(0,5);
    sigma_b12 ~ normal(0,5);
    sigma_d0 ~ normal(0,5);
    sigma_d1 ~ normal(0,5);
    sigma_rMid ~ normal(0, 2);
    sigma_rScale ~ normal(0, 2);
    
    //Random Effects
    b0_raw ~ normal(0, 1);
    b00_raw ~ normal(0, 1);
    b1_raw ~ normal(0, 1);
    b11_raw ~ normal(0, 1);
    b12_raw ~ normal(0, 1);
    d0_raw ~ normal(0, 1);
    d1_raw ~ normal(0, 1);
    rMid_raw ~ normal(0, 1);
    rScale_raw ~ normal(0,1);
}
