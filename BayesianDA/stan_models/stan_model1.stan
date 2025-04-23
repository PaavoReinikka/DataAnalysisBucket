// I'm writing STAN model from rstudio
      
      data {
        int N[4]; //sample sizes
        int y[4];  // deaths
        vector[4] x; // dozage
        vector[2] mu; // mu_alpha, mu_beta
        matrix[2,2] sigma; // variance-covariance matrix
      }
      
      parameters {
        vector[2] theta; // sampled alpha, beta
      }
      
      model {
        theta ~ multi_normal(mu, sigma);
        y ~ binomial_logit(N, theta[1] + theta[2]*x);
      }
      
      generated quantities {
      }
      
      
