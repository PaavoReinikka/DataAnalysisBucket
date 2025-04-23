# CALLs SOURCE FROM SA_functions
#devtools::install_github("b0rxa/scmamp")

require(ggplot2)
require(scmamp)

rank_diagram <- function(ranks, file_out = 'input', invert = T) {
  #' @title TODO: unfinished.
  r <- if(is.character(ranks)) {
    read.csv(ranks,row.names = 1)
  } else ranks
  r <- if(invert) t(r) else r
  if(file_out=='input') {
    file_out <- if(is.character(ranks)) paste(strsplit(ranks,'.csv')[[1]],'.jpeg',sep='') else 
      paste('from_console.jpeg', sep='')
    jpeg(file_out,width = 1400, height = 600)
    plotCD(r, cex = 0.4)
    dev.off()
  } else plotCD(r)
}

################### Function for ranking sobol indices. Same primary functionality as with coefficient rakings ###

rank_sobol <- function(ss, ties = 'max', coef_names = NULL, file_out = NULL) {
  #' @title Function for calculating ranks for sobol indices, such as they are when produced by SA_sobol.
  #' @description rank_sobol works in a similar way as rank_coefs, but needs SA_sobol$data as input.
  #' @param coef_names An additional parameter for choosing either subsets or changing the order of variables.
  #' This might be useful if SA_sobol is originally run with different ordering for plotting/subsetting purposes.
  #' @return Data frame of ranks, or outputs to separate files for 1st and total order inices.
  if(is.character(ss)) {
    file_in <- ss
    ss <- read.csv(ss)
  }
  if(is.null(coef_names)) coef_names <- names(ss)[-c(1,2)] 
  
  dss <- apply(ss[,coef_names] %>% abs(),1,
                function(data) {
                  r <- rank(-data, ties.method = ties)
                }) %>% as.data.frame()
  colnames(dss) <- paste(ss$Index, ss$Lambda)
  
  if(is.null(file_out)) dss else {
    if(file_out=='input') file_out <- strsplit(file_in,'.csv')[[1]]
    file_1 <- paste(file_out, '_S1_RANKS.csv', sep='')
    file_2 <- paste(file_out, '_ST_RANKS.csv', sep='')
    ind <- 1:ncol(dss)
    ds1 <- dss[,-which(ind%%2==0)]
    ds2 <- dss[,which(ind%%2==0)]
    write.csv(ds1, file = file_1, row.names = T)
    write.csv(ds2, file = file_2, row.names = T)
  }
}

rank_model <- function(mdl, ties = 'min', file_out = NULL, scale = T) {
  ret <- data.frame(coef(mdl)[-1,] %>% as.matrix() %>% abs() %>% rank(ties.method = ties))
  
  colnames(ret) <- paste('lambda_',mdl$lambda,sep = '')
  rownames(ret) <- rownames(coef(mdl)[-1,] %>% as.matrix())
  ret <- if(ties == 'max' & scale) ret - min(ret) else ret
  if(is.null(file_out)) ret else {
    file <- paste(file_out, '_model_RANKS.csv', sep='')
    write.csv(ret, file = file, row.names = T)
  }
}

################### Function for ranking model coefficient for different pairs (alpha, lambda) ###################
rank_coefs <- function(df_x, df_y, alphas = NULL, lambdas = NULL, ties = 'min', scale = T, file_out = NULL, folds = 10) {
  #' @title Function for calculating ranks for glmnet models.
  #' @description Ranks glmnet model coefficients. If ties == 'max' & scale == TRUE, the results are scaled so as to
  #' make zero coefficients have rank zero, and the most significant coefficient have the rank equaling the number of non-zero coeffs.
  #' @param df_x/df_y Clean data frames or matrices.
  #' @param alphas Values of elasticnet metaparam. If no values are provided, cv is used to find
  #' a value.
  #' @param lambdas Values for elasticnet param. for shrinkage. Defaults to cv {0, min, 1se}.
  #' @param ties How rank( ) function should handle ties.
  #' @param scale If ties == max, scale decides whether or not the results are scaled so that least values -> 0.
  #' @param file_out Name tag for .csv file where the results are saved.
  #' @return Data frame of ranks, or outputs to a file.
  
  x <- if(is.data.frame(df_x)) df_x %>% as.matrix() else df_x
  y <- if(is.data.frame(df_y)) df_y %>% as.matrix() else df_y
  
  alphas <- if(is.null(alphas)) {
    require(glmnetUtils)
    cat('Using cva.glmnet for finding alpha. Criteria used: min(cvm).')
    fit <- cva.glmnet(x, y)
    min <- Inf
    ind <- NULL
    for (i in 1:length(fit$modlist)) {
      new <- min(fit$modlist[[i]]$cvm)
      if(new<min) {
        min <- new
        ind <- i
      }
    }
    fit$alpha[ind]
  } else alphas
  
  res <- foreach(i=alphas, .combine = 'cbind') %do% {
    cv.fit <- cv.glmnet(x, y, alpha = i, nfolds = folds)
    l <- if(is.null(lambdas)) {
      c(0, cv.fit$lambda.min, cv.fit$lambda.1se)
    } else if(is.integer(lambdas) & length(lambdas)<=3) {
      c(0, cv.fit$lambda.min, cv.fit$lambda.1se)[lambdas]
    } else lambdas
    
    foreach(k=l, .combine = 'cbind') %do% {
      mdl <- glmnet(x,y,alpha=i, lambda=k)
      tmp <- coef(mdl)[-1,] %>% as.matrix() %>% abs()
      ret <- rank(-tmp, ties.method = ties) %>% as.data.frame()
      colnames(ret) <- paste('Lambda_',k,' Alpha_',i, sep = '')
      rownames(ret) <- rownames(coef(mdl)[-1,] %>% as.matrix())
      ret
    }
  }
  
  res <- if(scale & ties=='max') {
    temp <- rownames(res)
    res <- sapply(res[,1:length(res[1,])],function(data){
      (data - min(data))
    }) %>% as.data.frame()
    rownames(res) <- temp
    res
  } else res
  
  if(is.null(file_out)) res else {
    file <- file_out   # CHANGED 27.8   paste(file_out, '_C_RANKS.csv', sep='')
    write.csv(res, file = file, row.names = T)
  }
  
}
