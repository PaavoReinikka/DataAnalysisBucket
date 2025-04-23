
A0Lmin <- SA_sigma_normalized(train_x, train_y,lambda = as.integer(1),alpha = 0, LOO = T)
A0Lmin <- A0Lmin$coef %>% t() %>% as.data.frame()
A0Lmin <- cbind(Index='CF_norm', Lambda = 'min', A0Lmin)

A0L1se <- SA_sigma_normalized(train_x, train_y,lambda = as.integer(2),alpha = 0, LOO = T)
A0L1se <- A0L1se$coef %>% t() %>% as.data.frame()
A0L1se <- cbind(Index='CF_norm', Lambda = '1se', A0L1se)

A0 <- rbind(A0Lmin, A0L1se)

res <- SA_mara(train_x, train_y, lambdas = as.integer(2:3),alpha = 0,threads = 2, folds = length(train_x[,1]),sorted = F)
res$plot
res$data$Lambda <- c('min', '1se')

ALL <- rbind(A0, res$data)

SA_plot(ALL, file_out = 'A0_SIGMA_vs_S1.pdf')

ALLAs <- rbind(A0, A.25, A.5, A.75, A1)
final <- cbind(Alpha=c(0,0,.25,.25,.5,.5,.75,.75,1,1), ALLAs)

write.csv(final, file = 'SA_results/ALL_COEFF_NORMALIZED.csv', row.names = F)


##########################################################

x <- stan.dize(train_x)
x <- x$data
y <- stan.dize(train_y)
y <- y$data
rank_coefs(df_x = x, df_y = y,alphas = c(0,.25,.5,.75,1), lambdas = as.integer(2),ties = 'max', scale = F, file_out = "compiled_ranks/Lambda_min_NORM_COEF_RANKS.csv")
rank_coefs(df_x = x, df_y = y,alphas = c(0,.25,.5,.75,1), lambdas = as.integer(3),ties = 'max', scale = F, file_out = "compiled_ranks/Lambda_1se_NORM_COEF_RANKS.csv")
###########################################################


s1_min_ranks <- read.csv('SA_results/compiled_ranks/EXP0_2000_S1_RANKS.csv', row.names = 1)[-2]
s1_min_ranks <- cbind(s1_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.25_2000_S1_RANKS.csv', row.names = 1)[-2])
s1_min_ranks <- cbind(s1_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.5_2000_S1_RANKS.csv', row.names = 1)[-2])
s1_min_ranks <- cbind(s1_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.75_2000_S1_RANKS.csv', row.names = 1)[-2])
s1_min_ranks <- cbind(s1_min_ranks, read.csv('SA_results/compiled_ranks/EXP1_2000_S1_RANKS.csv', row.names = 1)[-2])
names(s1_min_ranks) <- c('a0','a025','a05','a075','a1')
write.csv(s1_min_ranks, file = 'SA_results/compiled_ranks/Lambda_min_S1_RANKS.csv')

st_min_ranks <- read.csv('SA_results/compiled_ranks/EXP0_2000_ST_RANKS.csv', row.names = 1)[-2]
st_min_ranks <- cbind(st_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.25_2000_ST_RANKS.csv', row.names = 1)[-2])
st_min_ranks <- cbind(st_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.5_2000_ST_RANKS.csv', row.names = 1)[-2])
st_min_ranks <- cbind(st_min_ranks, read.csv('SA_results/compiled_ranks/EXP0.75_2000_ST_RANKS.csv', row.names = 1)[-2])
st_min_ranks <- cbind(st_min_ranks, read.csv('SA_results/compiled_ranks/EXP1_2000_ST_RANKS.csv', row.names = 1)[-2])
names(st_min_ranks) <- c('a0','a025','a05','a075','a1')
write.csv(st_min_ranks, file = 'SA_results/compiled_ranks/Lambda_min_ST_RANKS.csv')

s1_1se_ranks <- read.csv('SA_results/compiled_ranks/EXP0_2000_S1_RANKS.csv', row.names = 1)[-1]
s1_1se_ranks <- cbind(s1_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.25_2000_S1_RANKS.csv', row.names = 1)[-1])
s1_1se_ranks <- cbind(s1_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.5_2000_S1_RANKS.csv', row.names = 1)[-1])
s1_1se_ranks <- cbind(s1_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.75_2000_S1_RANKS.csv', row.names = 1)[-1])
s1_1se_ranks <- cbind(s1_1se_ranks, read.csv('SA_results/compiled_ranks/EXP1_2000_S1_RANKS.csv', row.names = 1)[-1])
names(s1_1se_ranks) <- c('a0','a025','a05','a075','a1')
write.csv(s1_1se_ranks, file = 'SA_results/compiled_ranks/Lambda_1se_S1_RANKS.csv')

st_1se_ranks <- read.csv('SA_results/compiled_ranks/EXP0_2000_ST_RANKS.csv', row.names = 1)[-1]
st_1se_ranks <- cbind(st_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.25_2000_ST_RANKS.csv', row.names = 1)[-1])
st_1se_ranks <- cbind(st_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.5_2000_ST_RANKS.csv', row.names = 1)[-1])
st_1se_ranks <- cbind(st_1se_ranks, read.csv('SA_results/compiled_ranks/EXP0.75_2000_ST_RANKS.csv', row.names = 1)[-1])
st_1se_ranks <- cbind(st_1se_ranks, read.csv('SA_results/compiled_ranks/EXP1_2000_ST_RANKS.csv', row.names = 1)[-1])
names(st_1se_ranks) <- c('a0','a025','a05','a075','a1')
write.csv(st_1se_ranks, file = 'SA_results/compiled_ranks/Lambda_1se_ST_RANKS.csv')

