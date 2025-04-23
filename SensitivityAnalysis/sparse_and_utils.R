

#CALLS SOURCE FROM SA_functions

################### Functions to help calculate conditional probability for a selected sparse set of coefficients #####
## These only take into acount which predictors got selected and nothing else.

update_model_count <- function(d = NULL, mdl) {
  if(is.null(d)) d <- numvecdict()
  d$append_number(which(coef(mdl) != 0), 1)
  d
}

calculate_proportion <- function(d, mdl) {
  sum(d[[which(coef(mdl)!=0)]])
}

################### Function for ranking sobol indices. Same primary functionality as with coefficient rankings ###

################################################
#An util for normalized fitting, and inv. operarations

stan.dize <- function(df,center = T, scale = T, m_in=NULL, sd_in=NULL) {
  
  sd <- if(is.null(sd_in)){
    sapply(df,var) %>% sqrt() %>% t() %>% as.data.frame()
    row.names(sd) <- 'dev'
  } else sd_in
  
  m <- if(is.null(m_in)) {
    colMeans(df) %>% t() %>% as.data.frame()
    row.names(m) <- 'mean'
  } else m_in
  
  if(ncol(df)==1) {
    m <- colMeans(df)
    names(m) <- paste(names(df),'mean')
    sd <- sqrt(var(df[1:nrow(df),]))
    names(sd) <- paste(names(df),'dev')
  }
  
  ds <- scale(df, center = center, scale = scale) %>% as.data.frame()
  list(data=ds, m=m, sd=sd)
}

################################################
#alternative: coefplot, coefpath etc.

plot_coef <- function(fit, intercept = T, tolerance = 1e-5, sparse = T) {
  
  s_beta <- if(sparse) sparse_coef(fit = fit, intercept, tolerance)  %>% list() %>% as.data.frame() else {
    if(intercept) coef(fit) %>% as.matrix() %>% as.data.frame() else coef(fit)[-1,] %>% as.matrix() %>% as.data.frame()
  }
  names(s_beta) <- 'value'
  s_beta <- s_beta %>% mutate(variable=row.names(s_beta),value=value)
  
  ggplot(data=s_beta) + 
    geom_point(aes(x=variable,y=value)) + 
    theme(axis.text.x = element_text(angle = 90))
  
  
}

################################################

s_predict <- function(model,data,intercept=F, tolerance = 1e-5) {
  b <- sparse_coef(model,tolerance)
  sparse_predict(data,b,intercept)
}


##################################################

sparse_coef <- function(fit, intercept = T, tolerance = 1e-5) {
  
  b_full <- coef(fit)
  b_sparse <- b_full[which(abs(b_full)>tolerance)]
  names(b_sparse) <- rownames(b_full)[which(abs(b_full)>tolerance)]
  if (intercept) b_sparse else b_sparse[-1]
}

#################################################

sparse_predict <- function(data,coef,intercept = F) {
  nimet <- names(coef)[2:length(names(coef))]
  
  if (intercept) {
    ds_x <- data %>% select(nimet) %>% mutate("(Intercept)" = 1)
    ds_x <- ds_x[c("(Intercept)",nimet)]  
    ret <- (ds_x %>% as.matrix()) %*% (coef %>% as.matrix())
  } else {
    ds_x <- data %>% select(nimet)
    (ds_x %>% as.matrix()) %*% (coef[2:length(coef)] %>% as.matrix())
  }
  
}  
#######################################################

sparsify_input <- function(full_input, mdl) {
  full_input %>% select(names(sparse_coef(mdl))[-1])
}

##############  To remove sparse (~0) input columns  ############

remove_sparse_columns <- function(ds, tolerance = 0.99) {
  #removes columns where proportion of zeros is higher than tolerance
  which_sparse <- foreach(i=ds) %do% {
    mean(i==0)>0.99
  }
  
  names_sparse <- names(ds)[which_sparse %>% unlist]
  
  ds %>% select(-names_sparse)
}

##############################################

# These can be used with latinhypercube [0,1], to expand to [..range..] NOT AT USE AT THE MOMENT

#############utils#############################
ranges <- function(x) map(x,range)

scale_one <- function(r,data) {
  ret <- (r[2]-r[1])*data
  ret+r[1]
}

scaled <- function(ranges,data) {
  map2(ranges,data,scale_one)
}
###############################################

########scales LHS to input range's############
return_scaled <- function(ds,LHsample) {
  cpy <- rep(list(LHsample),length(ds[1,]))
  rng <- ranges(ds)
  scaled(rng,cpy)
}
###############################################

