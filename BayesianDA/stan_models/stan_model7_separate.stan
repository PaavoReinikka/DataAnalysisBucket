// I'm writing STAN model from rstudio
      
      data {
        int<lower=0> N;
        int<lower=0> K;
        int<lower=0, upper=K> x[N];
        vector[N] y;
      }
      parameters {
        vector[K] mu;
        vector<lower=0>[K] sigma;
      }
      transformed parameters {
        
      }
      model {
        //sigma ~ scaled_inv_chi_square(K-1,8);
        y ~ normal(mu[x], sigma[x]);
      }
      generated quantities {
        real ypred[6];
        ypred = normal_rng(mu,sigma);
        //real ypred_1;
        //real ypred_2;
        //real ypred_3;
        //real ypred_4;
        //real ypred_5;
        //real ypred_6;
        //ypred_1 = normal_rng(mu[1], sigma[1]);
        //ypred_2 = normal_rng(mu[2], sigma[2]);
        //ypred_3 = normal_rng(mu[3], sigma[3]);
        //ypred_4 = normal_rng(mu[4], sigma[4]);
        //ypred_5 = normal_rng(mu[5], sigma[5]);
        //ypred_6 = normal_rng(mu[6], sigma[6]);
      }
      
      
