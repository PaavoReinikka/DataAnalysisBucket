
library(rstan)
library(gdata)
library(bayesplot)
library(aaltobda)
data("bioassay")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

mu <- c(0,10)
s <- matrix(c(4,10,10,100),nrow = 2)

stan_data <- list(mu=mu, sigma=s, N=bioassay$n, x=bioassay$x, y=bioassay$y)


write("// I'm writing STAN model from rstudio
      
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
      
      transformed parameters {   // No need to use loop, if x is a vector
        real alpha = theta[1];
        real beta = theta[2];
      
      }
      
      model {
        theta ~ multi_normal(mu, sigma);
        y ~ binomial_logit(N, alpha + beta*x);
      }
      
      generated quantities {
      }
      
      ",
      
      "stan_model1.stan"
      
      )


stanc("stan_model1.stan")

stan_model1 <- "./stan_model1.stan"


fit <- stan(file = stan_model1, data = stan_data, warmup = 200, iter = 1200, chains = 6, thin = 1, init = 'random')
?stan


posterior <- extract(fit)
str(posterior)

plot(fit)
plot(posterior$alpha, posterior$beta)

warms_alpha <- posterior$alpha
warms_beta <- posterior$beta

Rhat(cbind(warms_alpha[1:1000],warms_alpha[1001:2000],warms_alpha[2001:3000],
            warms_alpha[3001:4000],warms_alpha[4001:5000],warms_alpha[5001:6000]))

Rhat(cbind(warms_beta[1:1000],warms_beta[1001:2000],warms_beta[2001:3000],
           warms_beta[3001:4000],warms_beta[4001:5000],warms_beta[5001:6000]))

ess_bulk(cbind(warms_alpha[1:1000],warms_alpha[1001:2000],warms_alpha[2001:3000],
           warms_alpha[3001:4000],warms_alpha[4001:5000],warms_alpha[5001:6000]))

ess_bulk(cbind(warms_beta[1:1000],warms_beta[1001:2000],warms_beta[2001:3000],
           warms_beta[3001:4000],warms_beta[4001:5000],warms_beta[5001:6000]))


shinystan::launch_shinystan(fit)


### Optional splitting to 12 mixed subsets
Rhat(split(warms_alpha,f=1:12) %>% unlist())
Rhat(split(warms_beta,f=1:12) %>% unlist())


ind <- which(posterior$theta[,2]>0)
ld50 <- -posterior$theta[ind,1]/posterior$theta[ind,2]

hist(ld50[which(ld50<0.4)],breaks = seq(from=-.8,to = .4,by = .01), xlab = 'LD50')
hist(ld50)



