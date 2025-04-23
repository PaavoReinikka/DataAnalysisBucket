library(tidyverse)
library(reshape2)
library(ggplot2)

s1_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_S1_RANKS.csv', row.names = 1) %>% t()
st_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_ST_RANKS.csv', row.names = 1) %>% t()
cf_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_COEF_RANKS.csv', row.names = 1) %>% t()

var_s1 <- apply(s1_l1se, 2, function(data) var(data))
var_s1 <- rank(-var_s1, ties.method  = 'max')
var_st <- apply(st_l1se, 2, function(data) var(data))
var_st <- rank(-var_st, ties.method = 'max')
var_cf <- apply(cf_l1se, 2, function(data) var(data))
var_cf <- rank(-var_cf, ties.method = 'max')

var_all <- rbind(var_s1,var_st,var_cf) %>% as.data.frame()
write.csv(var_all, 'SA_results/compiled_ranks/Lambda1se_COEF_S1_ST_var_RANK_over_alpha.csv', row.names = T)

res <- rank_diagram(file_out = '',r=var_all, invert = F)

sd_s1 <- apply(s1_l1se, 2, function(data) var(data) %>% sqrt())
sd_s1 <- rank(-sd_s1, ties.method  = 'max')
sd_st <- apply(st_l1se, 2, function(data) var(data)%>% sqrt())
sd_st <- rank(-sd_st, ties.method  = 'max')
sd_cf <- apply(cf_l1se, 2, function(data) var(data)%>% sqrt())
sd_cf <- rank(-sd_cf, ties.method  = 'max')

sd_all <- rbind(sd_s1,sd_st,sd_cf) %>% as.data.frame()
write.csv(sd_all, 'SA_results/compiled_ranks/Lambda1se_COEF_S1_ST_sd_RANK_over_alpha.csv', row.names = T)


rank_diagram('SA_results/compiled_ranks/Lambdamin_COEF_S1_ST_sd_RANK_over_alpha.csv', invert = F)
rank_diagram('SA_results/compiled_ranks/Lambda1se_COEF_S1_ST_sd_RANK_over_alpha.csv', invert = F)
rank_diagram('SA_results/compiled_ranks/Lambdamin_COEF_S1_ST_var_RANK_over_alpha.csv', invert = F)
rank_diagram('SA_results/compiled_ranks/Lambda1se_COEF_S1_ST_var_RANK_over_alpha.csv', invert = F)
######################################################################################

cf_lmin <- read.csv('SA_results/compiled_ranks/Lambda_min_COEF_RANKS.csv', row.names = 1) %>% t()
cf_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_COEF_RANKS.csv', row.names = 1) %>% t()

cf_lmin <- cbind(Alpha, Lambda=lmin, cf_lmin) %>% as.data.frame()
cf_l1se <- cbind(Alpha, Lambda=l1se, cf_l1se) %>% as.data.frame()

moe_min <- cf_lmin %>% melt(id.vars=c('Alpha','Lambda'))
moe_1se <- cf_l1se %>% melt(id.vars=c('Alpha','Lambda'))
moe <- rbind(moe_1se, moe_min)


Alpha <- c(0,.25,.5,.75,1)
lmin <- rep('min',5)
l1se <- rep('1se',5)

s1_lmin <- read.csv('SA_results/compiled_ranks/Lambda_min_S1_RANKS.csv', row.names = 1) %>% t()
s1_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_S1_RANKS.csv', row.names = 1) %>% t()

s1_lmin <- cbind(Alpha, Lambda=lmin, s1_lmin) %>% as.data.frame()
s1_l1se <- cbind(Alpha, Lambda=l1se, s1_l1se) %>% as.data.frame()

moe_min <- s1_lmin %>% melt(id.vars=c('Alpha','Lambda'))
moe_1se <- s1_l1se %>% melt(id.vars=c('Alpha','Lambda'))
moe <- rbind(moe_1se, moe_min)

st_lmin <- read.csv('SA_results/compiled_ranks/Lambda_min_ST_RANKS.csv', row.names = 1) %>% t()
st_l1se <- read.csv('SA_results/compiled_ranks/Lambda_1se_ST_RANKS.csv', row.names = 1) %>% t()

st_lmin <- cbind(Alpha, Lambda=lmin, st_lmin) %>% as.data.frame()
st_l1se <- cbind(Alpha, Lambda=l1se, st_l1se) %>% as.data.frame()

moe_min <- st_lmin %>% melt(id.vars=c('Alpha','Lambda'))
moe_1se <- st_l1se %>% melt(id.vars=c('Alpha','Lambda'))
moe <- rbind(moe_1se, moe_min)


p <- ggplot(moe, aes(Lambda))
p + geom_bar(aes(fill=Alpha, y=value), stat = 'identity') +
  facet_wrap(~variable)
  
######################################################################################

rank_difference_plot <- function(file_in, alphas = 1:5) {
  
  ds <- read.csv(file_in, header = T, row.names = 1)
  ds <- ds %>% t() %>% as.data.frame()
  ds <- cbind(Alpha = c(0,.25,.5,.75,1)[alphas], ds[alphas,])
  moe <- ds %>% melt(id.vars = 'Alpha')
  
  p <- ggplot(data = moe) + 
    geom_bar(aes(x=Alpha, weight=value), color='black', fill = 'blue') +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~variable)
  
  ggsave(filename = paste(strsplit(file_in,'.csv')[[1]],'.pdf', sep=''), plot = p, device='pdf',
         width = 30, height = 20)
  
}


#######################################################################################
rank_comparison_plot <- function(r1, r2, tag1, tag2, alphas = 1:5, file_out = 'NONAME.pdf', geom='point') {
  
  if(is.character(r1)) {
    r1 <- read.csv(r1, header = T, row.names = 1)
    r2 <- read.csv(r2, header = T, row.names = 1)
  }
  r1 <- r1 %>% t() %>% as.data.frame()
  r2 <- r2 %>% t() %>% as.data.frame()
  
  aapo <- c(0,.25,.5,.75,1)[alphas]
  type1 <- rep(tag1, length(alphas))
  type2 <- rep(tag2, length(alphas))
  
  r1 <- cbind(Alpha = aapo, type = type1,r1[alphas,])
  r2 <- cbind(Alpha = aapo, type = type2,r2[alphas,])
  
  if(geom=='point') {
    moe <- rbind(r1,r2) %>% melt(id.vars = c('Alpha', 'type'))
    
    p <- ggplot(data = moe) +
      geom_point(aes(x=Alpha,y=value, color=type)) + 
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~variable)
  } else {  
    moe_1 <- r1 %>% melt(id.vars = c('Alpha','type'))
    moe_2 <- r2 %>% melt(id.vars = c('Alpha','type'))
    
    p <- ggplot(mapping = aes(x=Alpha, y=value)) +
      geom_histogram(data=moe_1, stat = 'identity',aes(x=Alpha-0.05), fill='blue',alpha=0.3) +
      geom_histogram(data = moe_2,stat = 'identity', fill='red', alpha = 0.9) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~variable)
  }
  
  ggsave(file_out, plot = p, device = 'pdf',width = 30,height = 20)
}

##################################################################################



############## Sensitivity indices (plots) #######################################
source('SA_functions.R')
SA_plot('SA_results/0_Sobol_2000.csv')
SA_plot('SA_results/0.25_Sobol_2000.csv')
SA_plot('SA_results/0.5_Sobol_2000.csv')
SA_plot('SA_results/0.75_Sobol_2000.csv')
SA_plot('SA_results/1_Sobol_2000.csv')

