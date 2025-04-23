// I'm writing STAN model from rstudio
      
      data {
        int<lower=0> N;
        vector[N] y;
      }
      parameters {
        real mu;
        real<lower=0> sigma;
      }
      transformed parameters {
        
      }
      model {
        y ~ normal(mu, sigma);
      }
      generated quantities {
        real ypred;
        ypred = normal_rng(mu, sigma);
      }
      
      
