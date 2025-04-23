require(tidyverse,warn.conflicts = F,quietly = T)
require(glmnet,warn.conflicts = F,quietly = T)
require(pse,quietly = T)
require(sensitivity,quietly = T)
require(doParallel,quietly = T)
require(reshape2,quietly = T)
require(ggplot2,quietly = T)
require(docstring,quietly = T)
require(readr,quietly = T)
#require(randtoolbox,warn.conflicts = F,quietly = T)
#source('./sparse_and_utils.R')
#source('./rank_functions.R')
#source('./residual_sampling.R')

############################## SA UTILS ####################################################

qdata <- function(p, data) quantile(x=data, probs=p)


initRun <- function(mdl) {    ### can use sparse version also
  function(data) {
    predict(mdl,newx = data %>% as.matrix())
  }
}

#### sparse predict wrapper ####
initSparse <- function(mdl) {
  function(data) {
    s_predict(mdl,data, intercept = T,tolerance = 1e-5)
  }
}


parameter_sample <- function(df_x, sample_size = 1000, type = 'sobol', file_out = NULL ) {
  #' @title A function for producing SA samples.
  #' @param df_x A clean data frame from input space.
  #' @param type The type of method to be used in sampling. Options: c("sobol", "innergrid", "grid")
  seed <- list()
  for(i in names(df_x)) {
    temp <- range(df_x %>% select(i))
    eval(parse(text=paste('seed$',i,' <- c(',temp[1],',',temp[2],')')))
  }
  sobSet <- parameterSets(seed,samples = sample_size, method = type)
  sobDs <- sobSet %>% as.data.frame()
  names(sobDs) <- names(train_x)
  if(is.null(file_out)) sobDs else write.csv(sobDs, file = paste(file_out, '_type_', type, sample_size,'.csv'), row.names = F)
}

SA_sigma_normalized <- function(df_x,df_y, lambda, alpha = 1, ties = 'max', LOO = F) {
  sigma_y <- sqrt(var(df_y))
  sigma <- lapply(df_x, function(data) sqrt(var(data))/sigma_y)
  if(is.integer(lambda)) {
    cv.fit <- cv.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), alpha = alpha, nfolds = if(LOO) length(df_x[,1]) else 10)
    lambda <- c(cv.fit$lambda.min, cv.fit$lambda.1se)[lambda]
  }
  fit <- glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), alpha = alpha, lambda = lambda)
  res <- (coef(fit)[-1] %>% as.vector()) * (sigma %>% as.numeric())
  res <- res^2
  names(res) <- names(sigma)
  ranks = rank(-abs(res), ties.method = ties)
  p <- ggplot() + 
    geom_point(aes(x=names(res), y=res)) +
    theme(axis.text.x = element_text(angle = 90))
  list(coef = res, ranks = ranks, plot = p) 
}

############################## SA fast99 ##################################################

SA_fast <- function(df_x, df_y = NULL, file_out = NULL, lambdas = NULL, alpha = 1, sample_size = 2000, threads = 3,
                    family = 'gaussian', subSet = NULL, facet = TRUE) {
  #' @title Function for sensitivity indices and related plots.
  #' @description Produces S1, ST for possibly multiple penalty parameter lambda values, and 
  #' visualizes the results. If df_y = NULL only provides sampling of the input space (X1, X2). Similarly, if file = <filename>, 
  #' a <filename>_sample_size_....csv is written and nothing returned. NOTE: file output is single data frame with size==sample_size,
  #' where as output$X1 and output$X2 are both sample_size/2.
  #' @examples 
  #' simple_example()
  #' SA_fast(x, y, alpha=0.5, sample_size = 2000)
  #' @param df_x Clean data frame with predictors (no irrelevant columns).
  #' @param df_y Clean data frame for responses (...). Defaults to NULL.
  #' @param lambdas Values of lambda, for which glmnet model and SI's are calculated.
  #' @param alpha The alpha value to be used with cv.glmnet and glmnet.
  #' @param sample_size Size of the sample.
  #' @param threads Number of threads to be used in looping lambdas.
  #' @param subSet Subset of predictors to be plotted (returned SA data is still on all preddictors).
  #' @return Returns a data frame with SA data, and a plot. Form: list(data = fastList, plot = p).
  #' @author Paavo Reinikka, for thesis work.
  
  
  ## xlist is arg.list for SA foo's (since 'qdata' needs data for quantiles)
  xlist <- lapply(df_x,list)
  n <- names(xlist) 
  q <- rep('qdata',length(n))
  
  #if no df_y, only X1, X2 samples returned (to be used with other SA functions)
  if (is.null(df_y)) {
    f <- fast99(model = NULL, factors = n, sample_size,q=q, q.arg = xlist)
    if (is.null(file_out)) list(X1 = f$X[1:(sample_size/2),], X2 = f$X[(sample_size/2+1):sample_size,]) else
      write.csv(f$X, row.names = F, file = paste(file_out,'_fast_sample_', sample_size,'.csv', sep = ''))
  
  # Otherwise, S1, ST are plotted and returned
  } else {
    
    if (is.null(lambdas) | is.integer(lambdas)) {
      fit <- cv.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(),alpha=alpha, family = family)
      lambdas <- if(is.integer(lambdas)) c(0,fit$lambda.min,fit$lambda.1se)[lambdas] else c(0,fit$lambda.min,fit$lambda.1se)
    }
    
    x <- df_x %>% as.matrix()  # for more concise notation
    y <- df_y %>% as.matrix()

    lassoList <- foreach(lmd = lambdas) %do% {
      glmnet(x,y,alpha = alpha,lambda = lmd, family = family)
    }
    
    cl <- makeCluster(threads)
    clusterExport(cl, varlist=c("qdata", "initRun"), env=environment())
    registerDoParallel(cl)

    fastList <- foreach(lmd = lassoList, .packages = c('sensitivity','glmnet','magrittr'),.combine='rbind') %dopar% {
      
      #gives indexes,lambdas, id's
      f <- fast99(model=initRun(lmd),factors = n,sample_size,q=q,q.arg=xlist)#, M = 10)
      corner <- data.frame(Index=c('S1','ST'),Lambda=c(lmd$lambda,lmd$lambda))
      cbind(corner,rbind( if(is.na(f$D1/f$V)) 0 else f$D1/f$V,if(is.na(f$Dt/f$V)) 0 else f$Dt/f$V))
    }

    stopCluster(cl)

    # ALPHABETICAL ORDER
    colnames(fastList) <- c('Index','Lambda',n)
    fastList <-  fastList[,order(names(fastList))] 
    fastList <- fastList %>% select(c('Index','Lambda')) %>% cbind(.,fastList %>% select(-c('Index','Lambda')))
    
    
    #colnames(fastList) <- c('Index','Lambda',n)
    moe <- reshape2::melt(id.vars = c('Index','Lambda'),data = if(is.null(subSet)) fastList else fastList %>% select(c('Index','Lambda',subSet)))

    p <- ggplot(data = moe,mapping = aes(x=variable,y=as.numeric(value),color=Index)) + 
      geom_point(mapping = aes(x=variable,y=as.numeric(value),color=Index, shape = factor(Lambda))) + 
      theme(axis.text.x = element_text(angle = 90))
    p <- if(facet) p + facet_wrap(~Lambda) else p
  
    if (is.null(file_out)) list(data = fastList, plot = p) else
      write.csv(fastList, row.names = F, file = paste(file_out,'_fast_sensitivity_', sample_size,'.csv', sep = ''))

  }
}


############################ SA LHS (sparse) ################################################


SA_lhs <- function(df_x, model = NULL, size = 1000, nboot = 50, repetitions = 1) {
  #' @title Function produces pse::LHS object for sensitivity information
  #' @description Mostly for sparse EN-models. Can also be used without a model
  #' to get a latinhypercube sampling of the input space. NOTE: does not behave well
  #' if there are extremely sparse input columns (>0.95 proportion of zeros). NOTE: due to use of sparse_predict
  #' only family = 'gaussian' models produce correct results (will be fixed if need arise).
  #' @param df_x Input space data frame. Does not need to be clean/sparse.
  #' @param model 'elnet'/'glmnet' or other model with coef and predict methods.
  #' @param size Sampling size.
  #' @param nboot Size of confidence bootstrap.
  #' @param repetitions Repetitions on analysis for some visualizations (e.g., plotcv)
  #' @return Returns a LHS object, or when model = NULL a sampling to be used with tell etc. decoupled methods. 
  #' Not recommended to use as the only SA approach, rather as an additional way to observe specific dimensions.
  #' @author Paavo Reinikka, for thesis work.
  #' @examples 
  #' simple_example2()
  #' 
  #' res <- SA_lhs(train_x, final_model, size = 100)
  #' plotprcc(res)
  #' plotscatter(res)
  
  n <- if (is.null(model)) names(df_x) else names(sparse_coef(model)[-1])
  xlist <- lapply(df_x %>% select(n),list)
  q <- rep('qdata',length(n))
    
  LHS(model = if (is.null(model)) model else initSparse(model), factors = n, size,  q=q, q.arg = xlist, nboot = nboot, repetitions = repetitions)
  
}



######################### SA MARA (s1's only) ###########################################


SA_mara <- function(df_x=NULL, df_y, lambdas = NULL, alpha = 1, sample_size = 1000, threads = 3,
                    family = 'gaussian', subSet = NULL, sparse = F, file_in = NULL, file_out = NULL, sorted = T, folds = 10) {
  #' @title Function for 1st order sensitivity indices and related plots.
  #' @description Produces S1 for possibly multiple penalty parameter lambda values, and 
  #' visualizes the results. If df_y is 'glmnet','elnet','fishnet' model, function returns S1's for that (single) model.
  #' If df_x,df_y are data, they are used as inputs for glmnet. If file_in != NULL, it only means, that *sampling* is given,
  #' not that it will be used for fitting a model(s). Model **always** comes from data or is given as a ready glmnet/etc. model.
  #' Sparsity can be forced, but if sample is from sparse space, it has to be specifically -> sparse = TRUE.
  #' @examples 
  #' SA_mara(x, y, family = 'poisson')
  #' @param df_x Clean data frame with predictors (no irrelevant columns). If NULL, need 'file_in' to give </path/filename>
  #' to .csv file which holds lhs,fast99 or other sample for sensitivity calculation.
  #' @param df_y Clean data frame for responses. Or glmnet (...) model.
  #' @param lambdas Values of lambda, for which glmnet model and SI's are calculated.
  #' @param alpha The alpha value to be used with cv.glmnet and glmnet. Always needed unlike with SA_sobol and SA_fast(when sampling).
  #' @param sample_size Size of the sample.
  #' @param threads Number of threads to be used in looping lambdas. NOTE: shouldn't be more than length(lambdas).
  #' @param subSet Subset of predictors to be plotted (returned SA data is still on all preddictors).
  #' @param sparse Forced sparsity can be applied. Performance-wise, there is no need for sparsity.
  #' @param file_in/out Name/tag for .csv input/output. For input, dir path is needed.
  #' @param sorted Alphabetical ordering, defaults to TRUE.
  #' @return Returns a data frame with SA data, and a plot. Form: list(data = sobolListMara, plot = p).
  #' @author Paavo Reinikka, for thesis work.
  
  
  ## xlist is arg.list for SA foo's (since 'qdata' needs data for quantiles)
  if (!is.null(file_in)) {
  
    X <- read.csv(file_in)
    X <- if(sparse & !is.data.frame(df_y)) sparsify_input(mdl = df_y, X) else X
    
    n <- names(X)
    sample_size <- length(X[,1])/length(n)
  
  } else {
    if(is.null(df_x)) stop('If df_x == NULL, need file_in!')
    ds <- if(sparse & !is.data.frame(df_y)) sparsify_input(mdl = df_y, df_x) else df_x
    
    xlist <- lapply(ds, list)
    n <- names(xlist) 
    q <- rep('qdata',length(n))
  
    FAST <- fast99(model=NULL,factors = n, sample_size, q=q,q.arg=xlist)
    X <- FAST$X
    
  }
  if (sum(class(df_y) %in% c('elnet','glmnet','fishnet')) > 0) {
    
    s <- if (sparse) sobolmara(model=initSparse(df_y), X1 = X) else sobolmara(model=initRun(df_y), X1 = X)
    corner <- data.frame(Index='S1',Lambda=df_y$lambda)
    sobolListMara <- cbind(corner,t(s$S))
    
  } else {
    
    if (is.null(lambdas) | is.integer(lambdas)) {
      slice <- if(is.integer(lambdas)) lambdas else 1:3
      fit <- cv.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(),alpha=alpha, family = family, nfolds = folds)
      lambdas <- c(0,fit$lambda.min,fit$lambda.1se)[slice]
    }
    lassoList <- foreach(lmd = lambdas) %do% {
      glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), alpha = alpha, lambda = lmd,family = family)
    }
    
    cl <- makeCluster(threads)
    clusterExport(cl, varlist=c("qdata", "initRun"), env=environment())
    registerDoParallel(cl)
    
    sobolListMara <- foreach(lso = lassoList, .packages = c('sensitivity','glmnet','magrittr'), .combine = 'rbind') %dopar% {
      s <- sobolmara(model=initRun(lso), X1 = X)
      corner <- data.frame(Index='S1',Lambda=lso$lambda)
      cbind(corner,t(s$S))
    }
    stopCluster(cl)
  }
  
  # ALPHABETICAL ORDER
  colnames(sobolListMara) <- c('Index','Lambda',n)
  sobolListMara <-  sobolListMara[,if(sorted) order(names(sobolListMara)) else names(sobolListMara)]
  sobolListMara <- sobolListMara %>% select(c('Index','Lambda')) %>% cbind(.,sobolListMara %>% select(-c('Index','Lambda')))
  moeMara <- reshape2::melt(id.vars = c('Index','Lambda'),
                            data = if(is.null(subSet)) sobolListMara else sobolListMara %>% select(c('Index','Lambda',
                                                                                                     subSet[if(sorted) order(subSet) else subSet])))
  
  p <- ggplot(data = as.data.frame(moeMara)) + 
    geom_point(mapping = aes(x=variable,y=as.numeric(value),color=Index)) + 
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~Lambda) 
  
  
  if (is.null(file_out)) list(data = sobolListMara, plot = p) else write.csv(sobolListMara,
                                                                             paste(file_out,'_S1_mara_',sample_size,'.csv',sep = ''),row.names = F)
}


  
  ############################## Sobol2007 #######################################
  
  SA_sobol <- function(df_x = NULL, df_y, lambdas = NULL, alpha = NULL, sample_size = 1000, threads = 3, family = 'gaussian',
                       subSet = NULL, sparse = F, file_in = NULL, file_out = NULL, sorted = T, FF = NULL, folds = 10) {
    #' @title Function for 1st and total order sensitivity indices and related plots.
    #' @description Produces S1,ST for possibly multiple penalty parameter lambda values, and 
    #' visualizes the results. If df_y is 'glmnet','elnet','fishnet' model, function returns S1,ST's for that (single) model.
    #' @param df_x Clean data frame with predictors (no irrelevant columns).
    #' @param df_y Clean data frame for responses. Or glmnet (...) model.
    #' @param lambdas Values of lambda, for which glmnet model and SI's are calculated. If NULL, default = 0, min, 1se.,
    #' which are acquired either via cv.glmnet or, if alpha = NULL, cva.glmnet. Lambdas can also be integers 1:3 and are then used
    #' to slice cv/cva {0,min,1se}.
    #' @param alpha The alpha value to be used with cv.glmnet and glmnet.
    #' @param sample_size Size of the sample.
    #' @param threads Number of threads to be used in looping lambdas. NOTE: shouldn't be more than length(lambdas).
    #' @param subSet Subset of predictors to be plotted (returned SA data is still on all predictors).
    #' @param sparse Forced sparsity. SHould be used with large sample sizes.
    #' @param file_in/out Name/tag for .csv input/output. Input with dir path.
    #' @param sorted Whether or not the results are sorted alphabetically.
    #' @param FF For factor fixing. list(name=<colname>, value=<value set for column>). FF$value defaults to colMean.
    #' @return Returns a data frame with SA data, and a plot. Form: list(data = sobolList2007, plot = p).
    #' @author Paavo Reinikka, for thesis work.
    
    ##  I have made this pretty robust, so that it works in a reasonable way regardless the 
    ##  combination of inputs. For that reason I could not avoid multiple checks/conditions.
    ##  Apologies, I know it's torturous to try and read these...
    
    if(!is.null(FF) & is.data.frame(df_x)) {
      
      if(is.null(FF$name)) stop('FF$name need to be given for factor fixing')
      value <- if(is.null(FF$value)) colMeans(df_x %>% select(FF$name)) else FF$value
      cat(paste('Factor Fixing', FF$name))
      n <- names(df_x)
      df_x[,FF$name] <- value
    }
    
    lambdas <- if(is.null(lambdas)) as.integer(1:3) else lambdas
    if(sum(class(df_y) %in% c('elnet','glmnet','fishnet')) > 0) {}
    else if(is.integer(lambdas) | is.null(alpha)) {
      if(is.null(alpha)) {
        require(glmnetUtils)
        cat('Using cva.glmnet for finding alpha & lambda. Criteria used: min(cvm).')
        fit <- cva.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), nfolds = folds)
        min <- Inf
        ind <- NULL
        for (i in 1:length(fit$modlist)) {
          new <- min(fit$modlist[[i]]$cvm)
          if(new<min) {
            min <- new
            ind <- i
          }
        }
        alpha <- fit$alpha[ind]
        lambdas <- c(0,fit$modlist[[ind]]$lambda.min,fit$modlist[[ind]]$lambda.1se)[lambdas]
      } else {
        if(length(lambdas)>3) stop('Too many integer lambdas! Only 3 first are used')
        fit <- cv.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), alpha = alpha, nfolds = folds)
        lambdas <- c(0,fit$lambda.min, fit$lambda.1se)[lambdas]
      }
    }
    
    if (!is.null(file_in)) {
      
      X <- read.csv(file_in)
      X <- if(sparse & !is.data.frame(df_y)) sparsify_input(mdl = df_y, X) else X
      
      half <- sample(1:length(X[,1]),length(X[,1])/2, replace=F)
      n <- names(X)
      sample_size <- length(X[,1])/length(half)
      X1 <- X[half,]
      X2 <- X[-half,]
      
    } else {
      
      if(is.null(df_x)) stop('If df_x == NULL, need file_in!')
      ds <- if(sparse & !is.data.frame(df_y)) sparsify_input(mdl = df_y, df_x) else df_x
      
      ## xlist is arg.list for SA foo's (since 'qdata' needs data for quantiles)
      xlist <- lapply(ds,list)
      n <- names(xlist) 
      q <- rep('qdata',length(n))
      
      FAST <- fast99(model=NULL,factors = n, sample_size, q=q,q.arg=xlist)
      half <- sample(1:length(FAST$X[,1]),length(FAST$X[,1])/2, replace=F)
      X1 <- FAST$X[half,]
      X2 <- FAST$X[-half,]
    }
    
    if (sum(class(df_y) %in% c('elnet','glmnet','fishnet')) > 0) {
      
      s <- sobol2007(model=if (sparse) initSparse(df_y) else initRun(df_y), 
                     X1 = X1, X2 = X2)
      corner <- data.frame(Index=c('S1','ST'), Lambda=c(df_y$lambda,df_y$lambda))
      s <- cbind(corner,rbind(t(s$S),t(s$T)))
      sobolList2007 <- s
      
    } else {
      
      if(is.null(df_x)) stop('If df_x == NULL, need "elnet" model as df_y!')
      
      #if (is.null(lambdas)) {
      #  fit <- cv.glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(),alpha=alpha, family = family)
      #  lambdas <- c(0,fit$lambda.min,fit$lambda.1se)
      #}
      lassoList <- foreach(lmd = lambdas) %do% {
        glmnet(df_x %>% as.matrix(), df_y %>% as.matrix(), alpha = alpha,lambda = lmd,family = family)
      }
      
      cl <- makeCluster(if(threads>length(lambdas)) length(lambdas) else threads )
      clusterExport(cl, varlist=c("qdata", "initRun"), env=environment())
      registerDoParallel(cl)
      
      sobolList2007 <- foreach(lso = lassoList, .packages = c('sensitivity','glmnet','magrittr'), .combine = 'rbind') %dopar% {
        s <- sobol2007(model=initRun(lso), X1 = X1, X2 = X2)
        corner <- data.frame(Index=c('S1','ST'), Lambda=c(lso$lambda,lso$lambda))
        cbind(corner,rbind(t(s$S),t(s$T)))
      }
      stopCluster(cl)
    }
    
    #ALPHABETICAL ORDER OR NOT
    
    colnames(sobolList2007) <- c('Index','Lambda',n)
    sobolList2007 <-  sobolList2007[,if(sorted) order(names(sobolList2007)) else names(sobolList2007)] 
    sobolList2007 <- sobolList2007 %>% select(c('Index','Lambda')) %>% cbind(.,sobolList2007 %>% select(-c('Index','Lambda')))
    
    #colnames(sobolList2007) <- c('Index','Lambda',n)
    tempList <- if(is.null(subSet)) sobolList2007 else sobolList2007 %>% select(c('Index','Lambda',
                                                                                  subSet[if(sorted) order(subSet) else subSet]))
    moe2007 <- reshape2::melt(id.vars = c('Index','Lambda'),data = as.data.frame(tempList))
    
    
    p <-  ggplot(data = moe2007) + 
      geom_point(mapping = aes(x=variable,y=as.numeric(value),color=Index)) + 
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~Lambda)
    
    
    if (is.null(file_out)) list(data = sobolList2007, plot = p) else {
      data_file <- if(sparse) paste(file_out,'_Sobol_sparse_',sample_size,'.csv',sep = '') else paste(file_out,'_Sobol_', sample_size, '.csv',sep = '')
      plot_file <- if(sparse) paste(file_out,'_Sobol_sparse_',sample_size,'.pdf', sep = '') else paste(file_out,'_Sobol_', sample_size,'.pdf', sep = '')
      
      write.csv(sobolList2007, data_file, row.names = F)
      ggsave(plot_file, plot = p, device = 'pdf')
    }
  }
  
  
  SA_plot <- function(data_in, file_out = NULL, path = NULL, device = 'pdf', scale = 1, dpi = 300, width = 0.5, height = 20, units = 'cm') {
    #' @title Plotting utility for SA outputs.
    #' @description Plots SA results to file with reasonable defaults.
    #' @param data_in Can be either a filename, complete with path, or a data frame from SA_... function.
    #' @param file_out Filename for output, complete with path. As default uses input files naming (or 'from_console.pdf'
    #' if data_in was a data frame)
    #' @param path Path to file destination can be given explicitly (or implicitly in filename).
    #' @param device Filetype, 'pdf', 'jpeg', etc. what ever is supported on 'this' platform by ggplot2.
    #' @param scale/dpi/units Self-explanatory. For tweaking perhaps.
    #' @param width/height Only height should be adjusted unless for specific reasons. Width is calculated from height
    #' and the number of different lambdas (elnet penalty params.) in the SA data_in. In any case, width should be relative to 
    #' the number of subplots/facets in the plot (facet_wrap is  ~lambdas). 
    data <- if(is.character(data_in)) read.csv(data_in) else data_in
    
    moe <- reshape2::melt(id.vars = c('Index','Lambda'),data = data)
    
    
    p <- ggplot(data = moe) + 
      geom_point(mapping = aes(x=variable,y=as.numeric(value),color=Index)) + 
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~Lambda)
    
    width <- width*height*length(data[,1])
    
    if(is.null(file_out)) {
      file_out <- if(is.character(data_in)) paste(strsplit(data_in,'.csv')[[1]],'.', device,sep='') else paste('from_console.',device, sep='')
    }
    ggsave(file_out,plot = p,device = device, path=path, scale = scale, dpi = dpi, width = width, height = height, units = units)
    
  }
  
