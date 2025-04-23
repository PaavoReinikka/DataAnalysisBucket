#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


init_experiment <- function() {
  require(dplyr, quietly = T,warn.conflicts = F)
  ds <- read.csv('../data/dataset_y_oat_v1.csv')
  
  target <- 'production'
  idVars <- c('ELY','Year')
  omit <- names(ds)[c(65:68)]
  
  train_x <- ds  %>% select(-target,-idVars,-omit)   # training with everything here
  train_y <- ds  %>% select(target) %>% log()
  
  #EDIT ONLINE 28.8
  train_y <- train_y %>% scale() %>% as.data.frame()
  
  write.csv(train_x, file = './SA_clean_data/predictors.csv', row.names = F)
  write.csv(train_y, file = './SA_clean_data/responses.csv', row.names = F)
  cat('Clean predictors and responses are saved to directory /SA_main/SA_clean_data!\n')
}

init_sample <- function(tag, size) {
  source('./SA_functions.R')
  size <- as.numeric(size)
  file_in <- './SA_clean_data/predictors.csv'
  file_out <- paste('./SA_clean_data/',tag)
  x <- read.csv(file = file_in, header = T)
  sam <- SA_fast(df_x = x, sample_size = size, file_out = file_out)
}


#NOTE: changed 23.8
SA_indices <- function(alpha, size, method) {
  source('./SA_functions.R')
  source('./sparse_and_utils.R')
  source('./rank_functions.R')

  ds_x <- read.csv('./SA_clean_data/predictors.csv')
  ds_y <- read.csv('./SA_clean_data/responses.csv')
  alpha <- as.numeric(alpha)
  size <- as.numeric(size)
  folds <- length(ds_x[1,]) #for nfolds
  
  if(method == 'sobol') {
    cv.fit <- cv.glmnet(ds_x %>% as.matrix(), ds_y %>% as.matrix(), alpha = alpha, nfolds = folds)
    
    rank_coefs(ds_x,ds_y,alphas = alpha, lambdas = c(cv.fit$lambda.min, cv.fit$lambda.1se),
               file_out = paste('./SA_results/EXPonentiated/EXP',alpha,sep = ''))
    
    SA_sobol(df_x=ds_x,df_y=ds_y,lambdas = c(cv.fit$lambda.min, cv.fit$lambda.1se),alpha = alpha, sample_size = size, sparse = F,
           threads = 2, file_out = paste('./SA_results/EXPonentiated/EXP',alpha,sep = ''),sorted = F)
    
    ds <- read.csv(paste('./SA_results/EXPonentiated/EXP',alpha,'_Sobol_',size,'.csv',sep = ''))
    rank_sobol(ds,file_out = paste('./SA_results/EXPonentiated/EXP',alpha,'_',size,sep = ''))
  
  } else if(method=='mara') {
    cv.fit <- cv.glmnet(ds_x %>% as.matrix(), ds_y %>% as.matrix(), alpha = alpha, nfolds = folds)
    
    SA_mara(df_x=ds_x,df_y=ds_y,lambdas = c(cv.fit$lambda.min, cv.fit$lambda.1se),
            alpha = alpha, sample_size = size, sparse = F,
            threads = 2, file_out = paste('./SA_results/EXPonentiated/EXP',alpha,sep = ''),sorted = F)
  }
  
}

all_coef_ranks <- function() {
  source('./SA_functions.R')
  source('./rank_functions.R')
  
  ds_x <- read.csv('./SA_clean_data/predictors.csv')
  ds_y <- read.csv('./SA_clean_data/responses.csv')
  
  folds <- length(ds_x[1,])
  
  rank_coefs(ds_x,ds_y, alphas = c(0,0.25,0.5,0.75,1), file_out = './SA_results/ALLCOEFS', folds = folds)
}

SA_ols <- function(size) {
  source('./SA_functions.R')
  source('./sparse_and_utils.R')
  
  ds_x <- read.csv('./SA_clean_data/predictors.csv')
  ds_y <- read.csv('./SA_clean_data/responses.csv')
  size <- as.numeric(size)
  #folds <- length(ds_x[1,]) for nfolds
  
  SA_sobol(df_x=ds_x,df_y=ds_y,lambdas = 0,alpha = 1, sample_size = size, sparse = F, threads = 2, file_out = './SA_results/ols',sorted = F)
}


if(args[1]=='init') init_experiment()
if(args[1]=='sample') init_sample(args[2],args[3])
if(args[1]=='sobol') SA_indices(args[2], args[3], 'sobol')
if(args[1]=='ols') SA_ols(args[3])
if(args[1]=='mara') SA_indices(args[2], args[3], 'mara')
if(args[1]=='allcoefs') all_coef_ranks()

###############################################################################################################################################






wrap <- function() {
  
x <- read.csv('./SA_clean_data/predictors.csv',header = T, as.is = T)
y <- read.csv('./SA_clean_data/responses.csv',header = T, as.is = T)
#################################################################################################################################################
#################### FAST SAMPLE TO ./SA_clean_data/<FILE> (and indices to ./SA_results/<FILE2>) ################################################

SA_fast(df_x = x,file_out = './SA_clean_data/EXAMPLE',sample_size = 2000)
SA_fast(df_x = x, df_y = y, file_out = './SA_results/EXAMPLE',sample_size = 2000)

######################################### MARA TO FILE ##########################################################################################
SA_mara(df_x = x,df_y = y,sparse = T,alpha = 1,file_in = './SA_clean_data/EXAMPLE_fast_sample_2000.csv',file_out = './SA_results/S_EXAMPLE', sample_size = 300)

##################################################################################################################################################
################### SOBOL TO FILE. FROM FILE SAMPLES AS WELL AS FROM FIT'S #######################################################################

#since sobol has problem running full_sample, we can sparsify. 2 choises: Either we give it a elnet model and letit sparsify in function. Or
# we can sparsify (with that elnet model) the input space outside function (either sample or original data). 

### THERE ARE NO SETTING FILES, SO values of type s$alpha[1] won't work. They have to be put in manually here.
if(is.na(s$lambda[1])) {
  cv.fit <- cv.glmnet(x %>% as.matrix(), y %>% as.matrix(), alpha=s$alpha[1],nfolds = if(s$LOO[1]) length(x[,1]) else 10, grouped = if(s$LOO[1]) F else T)
  lambda <- cv.fit$lambda.min
  fit <- glmnet(x %>% as.matrix(), y %>% as.matrix(), alpha = s$alpha[1], lambda = lambda)
} else {
  fit <- glmnet(x %>% as.matrix(), y %>% as.matrix(), alpha = s$alpha[1], lambda = s$lambda[1])
}
##################################################################################################################################################

#option 1
# MODEL & SAMPLE FILE
SA_sobol(df_x = NULL ,df_y = fit ,alpha = s$alpha[1],threads = s$threads[1],
         file_in = './SA_clean_data/EXAMPLE_fast_sample_2000.csv',file_out = './SA_results/EXAMPLE',
         sparse = s$sparse[1], sample_size = 2000)

#option 2
sparse_x <- sparsify_input(full_input = x, mdl = fit)

# MODEL & SPARSE DATA
SA_sobol(df_x = sparse_x, df_y = fit ,alpha = s$alpha[1],file_out = './SA_results/EXAMPLE2_intercept',
         sparse = s$sparse[1], sample_size = 2000)

# Or we can run full, with just a single fit:
# MODEL & FULL DATA
SA_sobol(df_x = x, df_y = fit, alpha=s$alpha[1],file_out = './SA_results/EXAMPLE3',
         sparse = FALSE, sample_size = 2000)

# Or full path, with smaller sample_size:
# NO MODEL, JUST DATA (IF NO LAMBDAS ARE DESIGNATED, DEFAULT {0,min,1se})
SA_sobol(df_x=x, df_y=y, alpha=s$alpha[1],file_out = './SA_results/EXAMPLE4_path', sample_size = 500)

##################################################################################################################################################
##################### RESIDUALS TO FILE ##########################################################################################################

#RESIDUAL SAMPLING (FULL RESIDUALS +/- ADDED AS NOISE)
coef_bootstrap(x,y,strap_size = 3000, alpha=1, threads = 10, file = './SA_results/EXAMPLE5')
#PROPORTION OF RESIDUALS +/- ADDED AS NOISE
coef_bootstrap(x,y,strap_size = 3000, alpha=1, threads = 10, randomness = list(type='proportion',value=0.5),file = './SA_results/EXAMPLE6')
#PROPORTION OF SINGLE PREDICTOR ('area' in this case) RANDOMIZED
coef_bootstrap(x,y,strap_size = 3000, alpha=1, threads = 10, type = 'single',randomness = list(type='area',value=0.5),file = './SA_results/EXAMPLE7')

}
