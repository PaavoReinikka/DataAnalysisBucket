// I'm writing STAN model from rstudio
      
      data {
        int<lower=0> N; // number of data points
        vector[N] x; // observation year
        vector[N] y; // observation number of drowned
        real xpred; // prediction year
        real tau;
      }
      parameters {
        real alpha;
        real beta;
        real<lower=0> sigma; //upper bound changed to lower bound
      }
      transformed parameters {
        vector[N] mu;
        mu = alpha + beta*x;
      }
      model {
        beta ~ normal(0, tau); //Here's the added prior
        y ~ normal(mu, sigma);
      }
      generated quantities {
        real ypred;
        ypred = normal_rng(alpha+beta*xpred, sigma); //changed xpred
      }
      
      
