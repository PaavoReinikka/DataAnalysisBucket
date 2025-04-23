library(aaltobda)
library(tidyverse)

data("bioassay")


??dmvnorm


make_cov <- function(alpha_var, beta_var, corr) {
  tmp <- corr * sqrt(alpha_var * beta_var)
  matrix(c(alpha_var, tmp, tmp, beta_var), nrow = 2, byrow = T)
  
}

s <- make_cov(4,100, .5)
mu <- c(0,10)

alpha <- 3
beta <- 9

ss <- rmvnorm(1000,mu, s)
ss <- ss[which(ss[,2]>0),]

ld50 <- -ss[,1]/ss[,2]

p_log_prior <- function(alpha, beta) {
  m <- cbind(alpha, beta) %>% as.matrix(,nrow=length(alpha[,1]))
  dmvnorm(m, mu, s, log = T)
}

p_log_prior(alpha, beta)


p_log_posterior <- function(alpha, beta, x = bioassay$x, y = bioassay$y, n = bioassay$n) {
  
  p_log_prior(alpha, beta) + bioassaylp(alpha, beta, x,y,n)
  
}

p_log_prior(alpha, beta)
p_log_posterior(alpha, beta)


bioassay_posterior_density_plot(c(-4,4),c(-10,30), bioassay$x, bioassay$y, bioassay$n)
??bioassay_posterior_density_plot


alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)


log_importance_weights <- function(alpha, beta) {

  tmp <- p_log_posterior(alpha, beta) - p_log_prior(alpha, beta)
  #round(tmp, digits=2)
  tmp
}

normalized_importance_weights <- function(alpha, beta) {
  tmp <- (log_importance_weights(alpha, beta) %>% exp()) / sum(log_importance_weights(alpha, beta) %>% exp())
  #round(tmp, digits=3)
  tmp
}

s_eff <- function(alpha, beta) {
  w <- normalized_importance_weights(alpha, beta)
  1/sum(w^2)
  
}

posterior_mean <- function(alpha, beta) {
  
  w <- normalized_importance_weights(alpha, beta)
  tmp <- sum(w*alpha)
  c(tmp,sum(w*beta))
}


log_importance_weights(alpha, beta) %>% round(.,digits=2)
normalized_importance_weights(alpha, beta) %>% round(.,digits=3)

posterior_mean(alpha, beta)
s_eff(alpha, beta)

############################GRID

sam1 <- rmvnorm(5000, mu, s)
sam2 <- rmvnorm(5000, mu, s)
sam <- rbind(sam1, sam2)
sum(normalized_importance_weights(sam[,1],sam[,2])>0)

cA <- rep(A,each =length(B))
cB <- rep(B, length(A))

cW <- normalized_importance_weights(cA,cB)
sum(cW)
sample(1:length(cW), 1000, replace = F, prob = cW)
