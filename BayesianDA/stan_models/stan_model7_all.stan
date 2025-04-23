// I'm writing STAN model from rstudio
      
      data {
      //COMMON DATA
      
        int<lower=0> N;
        int<lower=0> K;
        int<lower=0, upper=K> x[N];
        vector[N] y;
      }
      parameters {
        //PARAMS FOR POOLED
        real mu_P;
        real<lower=0> sigma_P;
        
        //PARAMS FOR SEPARATE
        vector[K] mu_S;
        vector<lower=0>[K] sigma_S;
        
        //PARAMS FOR HIERACHICAL
        real mu_0_H;
        real<lower=0> sigma_0_H;
        vector[K] mu_H;
        real<lower=0> sigma_H;
      }
      transformed parameters {
        
      }
      model {
      //MODEL FOR POOLED (P)
        y ~ normal(mu_P, sigma_P); //uniform priors
        
      //MODEL FOR SEPARATE (S)
        sigma_S ~ scaled_inv_chi_square(K-1,8); //weak prior
        y ~ normal(mu_S[x], sigma_S[x]);
      
      //MODEL FOR HIERARCHICAL (H)
        mu_0_H ~ normal(90,50); // weak hyperprior
        sigma_0_H ~ scaled_inv_chi_square(K-1,8); // weak hyperprior
        sigma_H ~ cauchy(0,8); // weak prior
        mu_H ~ normal(mu_0_H,sigma_0_H); // prior
        y ~ normal(mu_H[x], sigma_H);
      }
      generated quantities {
        real ypred_H[6];
        real ypred_S[6];
        real ypred_P;
        
        //PREDICTION FOR POOLED
        ypred_P = normal_rng(mu_P, sigma_P);
        
        //PREDICTION FOR SEPARATE
        ypred_S = normal_rng(mu_S, sigma_S);
        
        //PREDICTION FOR HIERACHICAL
        ypred_H = normal_rng(mu_H,sigma_H);
      }
      
      
