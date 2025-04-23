
#CALLS SOURCE FROM SA_functions



require(docstring)
require(binaryLogic)
require(dict)

mock_entrophy <- function(ds_LOO) {
  #' @title A mock entrophy, assuming either totally uncorrelated predictors,
  #' or "ignorant observer". This function can not be used (alone at least) for sensitivity analysis' purposes.
  #' Rather, just to observe the behaviour or elnet from a selection point of view.
  ret <- ds_LOO[,-c(1,2,3)]
  ret <- ret[,which(ret!=0)]
  -sum(ret*log(ret))
}

single_LOO <- function(ds_x, ds_y, lambda, alpha = 1, output_type = 'raw') {
  #' @title An experimental function for probing lambda path in elnet models.
  #' @description This function executes a full LOO cross-validation for a single pair (<alpha>, <lambda>),
  #' and collects either I) 'raw' data from calculations or II) summary values. (see details below)
  #' @param ds_x/y A clean data frame with input/output.
  #' @param lambda/alpha Elnet parameters to be used for the cv.
  #' @param output_type Either 'raw' or anything else for summary.
  #' @return Raw output consists of all validation errors (as many as length(ds_y)), mean training errors
  #' relative to mean ols training error, and binary representation of the model at each cv step. Non-raw
  #' values are summaries: skewness of validation errors, mean of validation errors, variance of validation errors,
  #' mean of ols -relative errors (mean of rel. means -> no statistical significance, just of pure interest),
  #' and finally, the predictor probabilities (simple, unconditioned) i.e., mean of binary representations column-wise.
  x <- ds_x %>% as.matrix()
  y <- ds_y %>% as.matrix()
  corner <- data.frame(rep(0,2+length(names(ds_x)))) %>% t()
  
  for(i in 1:length(ds_x[,1])) {
    
    ols <- glmnet(x[-i,],y[-i,], alpha = 1, lambda = 0)
    fit <- glmnet(x[-i,],y[-i,], alpha = alpha, lambda = lambda)
    ols_training_MAD <- mean(abs(y[-i,] - predict(ols, x[-i,])))
    fit_training_MAD <- mean(abs(y[-i,] - predict(fit, x[-i,])))
    error <- y[i,] - predict(fit, x)[i]
    model.id <- coef(fit)[-1] %>% as.vector()!=0
    corner <- rbind(corner,c(error, fit_training_MAD/ols_training_MAD, as.binary(model.id,logic=T)  ))
  }
  corner <- corner[-1,] %>% as.data.frame()
  if(output_type=='raw') {
    colnames(corner) <- c('single_error','Ols.relative_training_MAD',names(ds_x))
    corner
  } else {
    colnames(corner) <- c('LOO_mean_error','Ols.relative_training_MAD',names(ds_x))
    ret <- colMeans(corner)  %>% t()
    ret <- cbind(LOO_error_skewness = e1071::skewness(corner[,1]), LOO_error_var = var(corner[,1]), ret) %>% as.data.frame()
    ret
  }
}


################### Functions to help calculate conditional probability for a selected sparse set of coefficients #####
## These only take into acount which predictors got selected and nothing else.

model_delta <- function(ds_x, ds_y, lambdas, alpha = 1) {
  #' @title A function for quantifying the change rate of the model (firts order difference).
  #' @description This function should be used in conjunction with probability paths.
  #' @param ds_x/y Clean data frames.
  #' @param lambdas A vector of lambdas. Should be the same as with probability paths.
  #' @return Returns a data frame with unscaled differences and corresponding lambdas.
  x <- ds_x %>% as.matrix()
  y <- ds_y %>% as.matrix()
  corner <- matrix( c(0,lambdas[1]), ncol = 2)
  for(i in 2:length(lambdas)) {
    fit1 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i-1])
    fit2 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i])
    delta <- (coef(fit2)[-1] %>% as.vector()!=0) - (coef(fit1)[-1] %>% as.vector()!=0)
    corner <- rbind(corner, c(sum(delta),lambdas[i]))
  }
  ret <- corner %>% as.data.frame()
  names(ret) <- c('Delta','Lambda')
  ret
}

model_sliding_delta <- function(ds_x, ds_y, lambdas, weights = c(1,2,2,1), alpha = 1) {
  #' @title A function for quantifying the change rate of the model (2nd order centered difference).
  #' @description This function should be used in conjunction with probability paths.
  #' @param ds_x/y Clean data frames.
  #' @param weights Weights to be used in sliding window (sliding averaging).
  #' @param lambdas A vector of lambdas. Should be the same as with probability paths.
  #' @return Returns a data frame with unscaled differences and corresponding lambdas.
  x <- ds_x %>% as.matrix()
  y <- ds_y %>% as.matrix()
  scale <- 1/sum(weights)
  corner <- matrix( c(0,lambdas[1], 0, lambdas[2]), ncol = 2, byrow = T)
  for(i in 3:(length(lambdas)-2)) {
    f_m1 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i-1])
    f_m2 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i-2])
    
    f <- glmnet(x,y,alpha=alpha, lambda=lambdas[i])
    
    f_p1 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i+1])
    f_p2 <- glmnet(x,y,alpha=alpha, lambda=lambdas[i+2])
    
    d_m1 <- (coef(f)[-1] %>% as.vector()!=0) - (coef(f_m1)[-1] %>% as.vector()!=0)
    d_m2 <- (coef(f)[-1] %>% as.vector()!=0) - (coef(f_m2)[-1] %>% as.vector()!=0)
    d_p1 <- (coef(f_p1)[-1] %>% as.vector()!=0) - (coef(f)[-1] %>% as.vector()!=0)
    d_p2 <- (coef(f_p2)[-1] %>% as.vector()!=0) - (coef(f)[-1] %>% as.vector()!=0)
    
    delta <- scale*(weights %*% c(sum(d_m2), sum(d_m1), sum(d_p1), sum(d_p2)))
    corner <- rbind(corner, c(delta,lambdas[i]))
  }
  corner <- rbind(corner, matrix( c(0,lambdas[length(lambdas)-1], 0, lambdas[length(lambdas)]), ncol = 2, byrow = T))
  ret <- corner %>% as.data.frame()
  names(ret) <- c('Delta.2nd','Lambda')
  ret
}



build_probability_path <- function(ds_x, ds_y, lambda, alpha = 1, sample_proportion = 1, strap_size = 1000,
                               replace = F, residual_data = NULL, output_type = 'raw', file_out = NULL, threads = 10) {
  #' @title A helper function to provide aid in stability selection and conditional probabilities for post-lasso inference.
  #' @description Bootstraps over lambda path, and gives probabilities for variables to be included in the model (or 'raw' 
  #' frequencies of appearances).
  #' @param ds_x/y Clean data frames.
  #' @param lambda A set of lambdas, a path over which the probabilities are calculated.
  #' @param sample_proportion The proportion of data to be used in bootstrap (if length(ds_y), the whole data is used -> less random)
  #' , or proportion of the original data to be perturbed with residuals (if length(ds_y), the whole data is perturbed -> more random)
  #' @param strap_size Number of runs for each lambda. Overall strap_size*length(lambda) runs.
  #' @param replace Self-evident. User should be cognizent how this effects the function behaviour and the nature of randomness.
  #' @param residual_data From e.g., ols fit. Are used for residual sampling on full data. Existence of these data dictates the
  #' meaning of sample_proportion and how it translates to randomness in bootstrap/res.sampling.
  #' @param threads Threads should not be more than length(lambda).
  #' @return Returns either a data frame or, if file_out given, outputs to file. Values are either raw frequencies or probabilities.
  
  cl <- makeCluster(threads)
  clusterExport(cl, varlist=c("build_binlist"), env=environment())
  registerDoParallel(cl)
  
  res <- foreach(i=1:length(lambda),.packages = c('glmnet','magrittr','binaryLogic','dplyr'), .combine = 'rbind') %dopar% {
    build_binlist(ds_x=ds_x, ds_y=ds_y, lambda= lambda[i], alpha=alpha, sample_proportion=sample_proportion,
                  strap_size=strap_size, replace=replace, residual_data = residual_data, output_type = output_type)
  }
  stopCluster(cl)
  if(!is.null(file_out)) {
    write.csv(res, file = paste(file_out, '.csv', sep = ''), row.names = F)
  } else res
}


build_binlist <- function(ds_x, ds_y, lambda, alpha = 1, sample_proportion = 1, strap_size = 1000, replace = F,
                          residual_data = NULL, output_type = 'raw') {
  #' @title Function for testing how many/which different models get selected over elnet.
  #' @description This is mostly called within build_probability_path, but can also be used for exploratory purposes
  #' and plotting the behaviour of variables for a single lambda within a bootstrap. Useful utility for plotting: 
  #' corrplot::corrplot( X, is.corr = F).
  #' @param args See build_probability_path( ). 
  
  x <- ds_x %>% as.matrix()
  y <- ds_y %>% as.matrix()
  len <- as.integer(sample_proportion*length(ds_x[,1]))
  tmp <- matrix(rep(0,length(ds_x[1,])),ncol = length(ds_x[1,]))
  
  if(is.null(residual_data)) {
    for(i in 1:strap_size) {
    
    #sample given proportion of training data, w/o replacement. 
      r <- sample(1:length(ds_x[,1]), len, replace = replace)
      x_r <- x[r,]
      y_r <- y[r,]
    # And fit (fit following same lambda as with the target i.e., "final" model.
      fit <- glmnet(x_r,y_r, alpha=alpha, lambda = lambda)
      tmp2 <- coef(fit)[-1] %>% as.vector()!=0
      
    # logic to binary to id value
      tmp <- rbind(tmp, as.binary(tmp2, logic = T))
    }
  } else {
    for(i in 1:strap_size) {
    #sample given proportion of residuals and add to responses.
      r <- sample(1:length(ds_x[,1]), len, replace = replace)
      signs <- sample(c(-1,1), length(r), replace=T)
      
      
      y_r <- ds_y + sample( c( rep(0,length(ds_x[,1])-len ), signs*residual_data[r,] ), replace = F)
    # And fit (fit following same lambda as with the target i.e., "final" model.
      fit <- glmnet(x, y_r %>% as.matrix(), alpha=alpha, lambda = lambda)
      tmp2 <- coef(fit)[-1] %>% as.vector()!=0
      
    # logic to binary to id value
      tmp <- rbind(tmp, as.binary(tmp2, logic = T))
    }
  }
  ret <- tmp[-1,] %>% as.data.frame()
  if(output_type=='raw') {
    names(ret) <- names(ds_x)
    ret <- ret %>% mutate(Lambda = lambda)
  } else {
    ret <- colMeans(ret)
    ret <- c(ret,lambda) %>% t() %>% as.data.frame()
    names(ret) <- c(names(ds_x),'Lambda')
  }
  ret
}


condition_groups <- function(ds_x, ds_y, lambda, alpha = 1, sample_proportion = 1, strap_size = 1000, replace = T) {
  #' @title This is an experimental function for calculating the frequencies for different combinations of variables 
  #' appearing together over elnet for a given alpha, lambda. 
  #' @description Bootstrap can be tuned to take subsamples w/o replacement, from the original data. The utility of this 
  #' function is mostly to give manual tool for probing the appearance of a given model in a bootstrap. Functions
  #' build_binlist and it's parallel "expansion" build_probability_path are much more versatile functions for
  #' model selection/validation.
  #' @param args All arguments are obvious, excluding maybe sample_proportion, which just defines the proportion of data
  #' to be used in a single fit within the bootstrap.
  #' @return Returns a numvecdict, where the keys are variable indeces (according to input design data's naming order)
  #' and values are vectors of 1's (for each appearance in a fit within bootstrap). Specific probabilities can be calculated 
  #' with calculate_proportion.

  corner <- numvecdict()
  
  for(i in 1:strap_size) {
    #sample given proportion of training data, w/o replacement. CV best min/1se lambda
    # (assumed one of those being the method used for fitting/selecting the model for which we want some conditioning)
    
    r <- sample(1:nrow(ds_x), as.integer(sample_proportion*nrow(ds_x)), replace = replace)
    
    #choose either lambda.min:= 1 or lambda.1se:= 2
    # And fit (fit following same heuristic as with the target i.e., "final" model.
    fit <- glmnet(ds_x[r,] %>% as.matrix(), ds_y[r,] %>% as.matrix(), alpha = alpha, lambda = lambda)
    
    #update dictionary with "selected" model.
    corner$append_number(which(coef(fit)!=0), 1)
  }
  #append the strap_size for proportion calc.
  corner$append_number(c(0,0,0), strap_size)
  corner
}
docstring(condition_groups)

calculate_proportion <- function(d, target_mdl) {
  sum(d[[which(coef(target_mdl)!=0) %>% as.numeric()]]) / d[[c(0,0,0)]]
}

############################################################################################
############################################################################################
###########################################################################################

fix_confidence <- function(predictions, p = c(0.05,0.95)) {
  #' A function to fix predictions from uncertainty predict. It gives p quantiles of predictions.
  #' The returned value is easy to e.g., cbind to an existing df with possible actual predictions and id tags.
  #' @param predictions Is a matrix from uncertainty_predict, or a data frame. 
  #' p Is at least of length 2 el. 0...1
  #' @return Returns <p> named columns for each prediction point.
  q_y <- apply(predictions[,1:length(predictions[1,])], 2, FUN = function(data) {
    quantile(data, probs = p)
  })
  
  t(q_y)
  
}
#EXAMPLE:

#res <- coef_bootstrap(train_x,train_y,threads = 10,strap_size = 10,type = 'permute', lambda=0.002)

#pred_y <- uncertainty_predict(test_x[1:length(test_x[,1]),],res$lasso)



uncertainty_predict <- function(input, coef_ds, intercept = T) {
  #' @title A function for point-wise uncertainty estimates
  #'
  #' @param Input Is clean (no excess columns) input data frame at the point(s) at which we wish to predict the output.
  #' @param coef_ds Comes from some parametric estimation process of the model coefficients, or from bootstrapping.
  #' @param intercept Whether or not model has an intercept or not, depends on coef_bootstrap(..., intercept = T). 
  #' @return Fuction calculates as many estimates as there are samples of coefficients in coef_ds (rows).
  
  x <- if(intercept) cbind("(Intercept)" = 1,input) else input
  (coef_ds %>% as.matrix()) %*% (x %>% as.matrix() %>% t())
}

########################################################################################################
########   this is for randomizing a column at a time (can add other types if necessary)  ##############
########################################################################################################

randomize_column <- function(ds, cname, randomness = 1.00, type = 'range', argList = NULL) {
  #' @title A function for randomizing a single column (by "name") in a data frame.
  #' @description Mostly an utility function to help observe the effects of a single predictor on an output.
  #' @param ds Data frame. 
  #' @param cname Name of the column to be randomized.
  #' @param randomness The proportion of random elements in the randomized column (default 1.0)
  #' @param type Types "range" (default) and "rnorm" .
  #' @param argList Supported with defaults, otherwise args = list(min,max) / list(mean,sd)
  #' @return A new data frame with column "cname" randomized.
  
  l <- length(ds[,cname])
  n <- as.integer(randomness*l)
  ind <- sample(1:l,n,replace=FALSE)
  ds[,cname] -> newcol
  
  if (type=='range') {
    
    newcol <- if (is.null(argList)) runif(l,min=min(ds[,cname]), max=max(ds[,cname])) else runif(l,min=argList$min, max=argList$max)
    ds[ind, cname] <- newcol[ind]
    
  } else if (type == 'rnorm') {
    
    newcol <- if (is.null(argList)) rnorm(l,mean=mean(ds[,cname]), sd = sd(ds[,cname])) else rnorm(l,mean=argList$mean, sd = argList$sd)
    ds[ind, cname] <- newcol
    
  }
  ds
}


########################################################################################################
####### This produces the bootstrapped coefs ###########################################################
########################################################################################################

coef_bootstrap <- function(x, y, strap_size = 1000, alpha = 1, lambda = NULL, family = 'gaussian', type = 'residual',
                           randomness = list(type='amplitude', value=1.0), threads = 4, subSet = NULL, intercept = T,
                           file = NULL) {
  #' @title Function to provide bootstrapped coefficients and spread / variance plots for them. 
  #' @description 
  #' If coefficient variances are needed, they should be calculated using the coefficient samples that this function returns.
  #' This function only takes clean data frames as x, y.
  #' NOTE: By giving a predictor by name to parameter 'randomness$type', and type = 'single'
  #' only that one predictor is randomized. This can be used to
  #' probe sensitivity/incertainty around an individual observation, caused by single noisy predictor. 
  #' @param x The predictor data frame
  #' @param y The response data frame
  #' @param strap_size Defaults to 1000, means simply the size of bootstrap...
  #' @param lambda Value to be used for glmnet. NOTE: should be left NULL (default uses cv min)
  #' @param type Either 'residual' for residual sampling, or 'permute' for response permutations. 'single' for only a single predictor. 
  #' @param randomness List(type,value), where type can be 'amplitude' (for residual noise) or 'proportion' (for proportion of
  #' responses that get added noise). For proportion, value must be in 0...1. NOTE 2: 
  #' randomness amount can only be defined for residual, or single column versions -> for single column version 
  #' randomness$type= <colname> and randomness$value is the amount of randomness (0...1). 
  #' @param subSet Chooses a specific subset of predictors to be plotted (to aid e.g., FA)
  #' @param intercept Whether or not include intercept in the model. Uncertainty_predict should have corresponding intercept flag.
  #' @param file A filename where results are saved. If given, nothings is returned.
  #' @return list(ols=ols_coefficients, lasso=lasso_coefficients, spread=plot). Nothing, if filename is given.
  #' If file = <filename>, outputs <filename>_<type>_..._.csv file, where information about type and amount of randomness if enbedded in filename.
  #' @examples 
  #' res <- coef_bootstrap(train_x, train_y)
  #' res$spread
  
  base_mdl <- glmnet(x %>% as.matrix(), y %>% as.matrix(), alpha = 1, lambda = 0, family = family)
  base_y <- predict(base_mdl, newx = x %>% as.matrix(), type = 'response')
  predictors <- names(x)
  
  res <- abs(y-base_y)
  l <- length(res[,1])
  
  a <- if(randomness$type =='amplitude') randomness$value else if (randomness$type=='proportion') {
    rand_amount <- as.integer(l*randomness$value)
    sample(c(rep(1,rand_amount),rep(0,l-rand_amount)),l,replace=F)
  } else 1.0
  
#  registerDoParallel(threads)
  
  cl <- makeCluster(threads)
  if (type=='single') clusterExport(cl, varlist=c("randomize_column"), env=environment())
  registerDoParallel(cl)
  
  if (type =='residual') {
    cat('residual...')
    residual_coefs <- foreach(i=1:strap_size, .combine = 'cbind', .packages = c('glmnet','magrittr')) %dopar% {
        
      # sample +/- from residuals and add to y
      t <- sample(c(-1,1), length(res[,1]), replace=T)
      perturb <- t*res[sample(1:l,l,replace = F),]
      
      # If randomness[1]=='amplitude', multiply with that (noise amplitude). Otherwise
      # multiply 'randomness[2]' proportion of these with ones, rest with zeros (in random order)
      y_rand <- a*perturb + y
      
      # This is not ideal, but can be fixed later (now it just slows down the loop)
      #if(family=='poisson') y_rand[which(y_rand<0),] <- -y_rand[which(y_rand<0),]
      
      fit <- cv.glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=alpha, family = family)
      olsModel <- glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=1, lambda = 0, family = family)
      lassoModel <- glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=alpha, 
                         lambda = if (is.null(lambda)) fit$lambda.min else lambda, family = family)
      
      rbind(coef(olsModel) ,coef(lassoModel))
    }
    cat('Done!')
  } else if (type == 'permute') {
    cat('permute...')
    residual_coefs <- foreach(i=1:strap_size, .combine = 'cbind', .packages = c('glmnet','magrittr')) %dopar% {
        
      y_rand <- y[sample(1:length(res[,1]), replace = F),]
        
      fit <- cv.glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=alpha, family = family)
      olsModel <- glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=1, lambda = 0, family = family)
      lassoModel <- glmnet(x %>% as.matrix(), y_rand %>% as.matrix(),alpha=alpha,
                           lambda = if (is.null(lambda)) fit$lambda.min else lambda, family = family)
        
      rbind(coef(olsModel) ,coef(lassoModel))
    }
    cat('Done!')
  } else if (type == 'single') {
    cat('single...')
    residual_coefs <- foreach(i=1:strap_size, .combine = 'cbind', .packages = c('glmnet','magrittr')) %dopar% {
      
      x_rand <- randomize_column(x, randomness$type, randomness$value)
      
      fit <- cv.glmnet(x_rand %>% as.matrix(), y %>% as.matrix(),alpha=alpha, family = family)
      olsModel <- glmnet(x_rand %>% as.matrix(), y %>% as.matrix(),alpha=1, lambda = 0, family = family)
      lassoModel <- glmnet(x_rand %>% as.matrix(), y %>% as.matrix(),alpha=alpha,
                           lambda = if (is.null(lambda)) fit$lambda.min else lambda, family = family)
      
      rbind(coef(olsModel) ,coef(lassoModel))
    }
    cat('Done!')
    
  } else stop('type arg. has to be one of "permute", "residual" or "single"!')
  
  stopCluster(cl)
    
  half <- 1:(1+length(predictors))
  alls <- residual_coefs %>% as.matrix()
  if (is.null(subSet)) subSet <- (if (intercept) c('(Intercept)',predictors) else predictors)
  
  
  olsBeta <- t(alls[half,]) %>% as.data.frame() %>% mutate(model=rep('ols',length(.[,1]))) %>% select('model',subSet)
  lassoBeta <- t(alls[-half,]) %>% as.data.frame() %>% mutate(model=rep('lasso',length(.[,1]))) %>% select('model',subSet)
  
  ##########coefficient samples to be returned############
  ols_coef_active <- t(residual_coefs[half,]) %>% as.matrix() %>% as.data.frame()
  lasso_coef_active <- t(residual_coefs[-half,]) %>% as.matrix() %>% as.data.frame()
  row.names(lasso_coef_active) <- 1:strap_size
  row.names(ols_coef_active) <- 1:strap_size
  ########################################################
  
  #####coefficient quantiles and means to be returned#####
  ols_range <- sapply(X = ols_coef_active, function(data) {
    
    q <- quantile(data, probs = c(0.05, 0.95))
    rbind(q[1],if (q[1]==0 & q[2]==0) 0 else mean(data),q[2])
  
  }) %>% as.data.frame()
                      
  lasso_range <- sapply(X = lasso_coef_active, function(data) {
    
    q <- quantile(data, probs = c(0.05, 0.95))
    rbind(q[1],if (q[1]==0 & q[2]==0) 0 else mean(data),q[2])
    
  }) %>% as.data.frame()
  
  rownames(ols_range)<- c('min','mean','max')
  rownames(lasso_range)<- c('min','mean','max')
  #######################################################
  
  ###########molten data for ggplot###############
  moe_ols <- olsBeta %>% melt(id.vars = 'model')
  moe_lasso <- lassoBeta %>% melt(id.vars = 'model')
  moe <- rbind(moe_ols,moe_lasso)
  
  moe2 <- ols_range %>% mutate(measure=rownames(.)) %>% melt(id.vars = 'measure')
  moe3 <- lasso_range %>% mutate(measure=rownames(.)) %>% melt(id.vars = 'measure')
  
  ##################plot#########################
  p1 <- ggplot(data = moe) +
    geom_boxplot(mapping = aes(x=variable,y=value,color=model)) + 
    ggtitle(paste('Coefficient spread over bootstrap of size ',strap_size)) +
    theme(axis.text.x = element_text(angle = 90))
  
  
  p2 <- ggplot(data=moe2) + 
    geom_point(aes(x=variable,y=value,color=measure)) +
    ggtitle(paste('Ols coefficient 0.05/mean/0.95 over bootstrap of size ',strap_size)) +
    theme(axis.text.x = element_text(angle = 90))
  
  p3 <- ggplot(data=moe3) + 
    geom_point(aes(x=variable,y=value,color=measure)) +
    ggtitle(paste('Elnet coefficient 0.05/mean/0.95 over bootstrap of size ',strap_size)) +
    theme(axis.text.x = element_text(angle = 90))
  
  ############# RETURN DATA AND PLOT, OR WRITE DATA TO FILE ##############################
  
  if (is.null(file)) list(ols=ols_coef_active, lasso=lasso_coef_active, ols_q = ols_range,
                          lasso_q=lasso_range, spread=p1, ols_plot = p2, lasso_plot = p3) else {
    
    #FILE NAMING
    file <-  if(type=='permute') paste(file, '_',type,'.csv',sep = '') else paste(file, '_',type,'_',randomness$type,'_',randomness$value,'.csv',sep = '')
    file2 <- paste('Q_ols_',file,sep = '')
    file3 <- paste('Q_elnet_',file,sep = '')
    
    write.csv(file = file, rbind(olsBeta, lassoBeta), row.names = F)
    write.csv(file = file2, ols_range,row.names = F)
    write.csv(file = file3, lasso_range,row.names = F)
  }
}


#############################################################################################################
#############################################################################################################

