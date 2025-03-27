picdgi_atitey <- function(dataset){
  
  library(HDInterval)
  
  ################# resampling function ---------------
  ################# 
  
  resample <- function(x, size, replace = TRUE, prob = NULL)
  {
    if(length(x)<1)
      if(!missing(size) && size>0)
        stop("Requested sample of size ", size, " from list of length 0")
    else
      x[FALSE]
    else if(length(x)==1)
    {
      if(missing(size) || size==1)
        x
      else if(size>=1 && replace==TRUE)
        rep(x, size)
      else if(size < 1)
        x[FALSE]
      else
        stop("Cannot cannot take a sample larger than the population",
             " when 'replace = FALSE'")
    }
    else
      sample(x, size, replace, prob)
  }
  
  ################# Computation of the AR and MA orders ---------------
  ################# 
  arma.data <- dataset[1,2:dim(dataset)[2]]
  arma.data <- as.numeric(arma.data)*10
  random_sample <- runif(17, min = min(as.numeric(arma.data)), max = max(as.numeric(arma.data)))
  arma.data <- as.numeric(c(arma.data,random_sample))
  coef <- 1.15 #1.12 ## 1/(ar.order+ma.order) #
  
  ################# Computation of the correlation matrix using the Browmian Law ---------------
  ################# 
  N <- 1000
  n <- dim(data.frame(dataset))[2]
  H = 0.9
  
  nu.square = 1.78
  mu = rep(0,N)
  mu <- data.frame(mu)
  mu <- t(t(mu))
  par=c(H,nu.square,coef)
  
  cor <- sapply(1:N, function(j) 0.5*(((abs(j+1))^(2.*par[1]))-2*((abs(j))^(2.*par[1]))+((abs(j-1))^(2.*par[1]))))
  cor.matrix <- toeplitz(cor)
  sigma <- par[2]*cor.matrix # variance of the mulvariate Normal distribution
  noise.dist <- rnorm(mu,sigma) # prior non stationnary states
  
  ################# arma model ---------------
  #################
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% arma coefficient
  #dataset <- 100*as.numeric(dataset[2:length(dataset)])
  dataset <- arma.data
  
  obs.vec <- sapply(1:n, function(i) dataset[i]) # observation vector
  obs.mat <- sapply(1:n, function(i) obs.vec[i]*t(obs.vec[i]))
  gam.fun <- sapply(1:n, function(i) abs((1/par[3])*obs.mat[i])) # coef.^-1.*ObsMat(1,j); % R^-1(y*y')
  gam.fun[is.na(gam.fun)] <- 0
  arma.coef.var <- sapply(1:n, function(i) exp(sum(diag(gam.fun[i]))))  # prior distribution P = exp(Tr(R^-1(y*y')))
  
  for(i in 1:n){
    arma.coef <- exp(sum(diag(gam.fun[i]))) # initial value of the distribution of ARMA coef P = exp(Tr(R^-1(y*y')))
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% smc
  sir.est <- c()
  variance <- matrix(0,nrow=n,ncol=N)
  inv.cov.mat <- matrix(0,nrow=n,ncol=N)
  particles = matrix(0,nrow=n,ncol=N)
  
  ## Calculate the prior density %%%%%%%%%%%%%%%
  for(t in 1:n){
    for(k in 1:N){
      particles[,k] <- arma.coef*noise.dist[k] # importance density is the transition density 
    }
    
    ## Calculate importance weights %%%%%%%%%%%%%%% 
    obs.vec[t] <- dataset[t] 
    variance[t,] = exp(particles[t,]) # exp(P53(1,t)./2)*vEps(t,1); % observation vector
    inv.cov.mat[t,] = 1/abs(variance[t,]) # inverse of the covariance matrix
    
    weights.pred = rep(1,N)/N
    weights.pred <- data.frame(weights.pred)
    weights_pred <- t(t(weights.pred))
    log.weight.0 = log(weights.pred)
    log.weight.1 = 0.0015; 
    log.weight.2 = -0.5*n*log(2*pi);
    log.weight.3= -0.5*log(abs(variance[t,]))
    log.weight.4 = -0.5*(obs.vec[t])*inv.cov.mat[t,]*obs.vec[t]
    log.weight.5 = sum(t(obs.vec[t] - particles[t,])*arma.coef.var[t]*(obs.vec[t] - particles[t,]))
    log.weight = log.weight.0 + log.weight.1 + log.weight.2 + log.weight.3 + log.weight.4 + log.weight.5
    
    ## Stabilize the importance weight %%%%%%%%%%%%%%% 
    d.max.weight <- max(log.weight)
    v.weight <- exp(log.weight - d.max.weight)
    v.weight <- v.weight/sum(v.weight)
    #v.weight.distribution <- v.weight
    sir.est[t] = sum(v.weight*particles[t,]);
    
    ## Compute the ESS %%%%%%%%%%%%%%%
    n.thr = 0.25*N
    n.eff = 1/(sum(v.weight^2))
    
    ## Resample and reset the weights %%%%%%%%%%%%%%% 
    if(isTRUE(n.eff < n.thr)==TRUE){ 
      index <- resample(v.weight,1:N)  
      index <- as.numeric(index$weights.pred)
      #particles <- particles[,index]
      particles <- matrix(index,             # Duplicate vector in matrix rows
                          nrow = n,
                          ncol = length(index),
                          byrow = TRUE)
      v.weight = rep(1,N)/N
      #v.weight.distribution <- v.weight
    } 
  }  
  
  ################# remove negative values ---------------
  #################  
  for(l in 1:n){
    if(isTRUE(sir.est[l] < 0)==TRUE){ 
      #if(sir.est[l] < 0) {
      sir.est[l] <- 0
    } else {
      sir.est[l] <- sir.est[l]
    }
    
  }
  
  sir.est[is.na(sir.est)] <- 0
  
  ################# IQR ---------------
  ################# 
  density.function <- as.numeric(v.weight$weights.pred)
  density.function[is.na(density.function)] <- 0
  density.new <- density.function[ !density.function <= exp(-20)] # remove values which are less than exp(-20) 
  score <- hdi(density.new, prob=95)
  
  ################# score ---------------
  ################# 
  score.min <- score[1]
  score.max <- score[2]
  hdi.range <- density.new[density.new >= score.min & density.new <= score.max]
  hdi.score <- (hdi.range/sd(hdi.range))^2
  sco.HDI <- mean(hdi.score)
  
  return(list(sir.est=sir.est, obs.vec=obs.vec, hdi.range=hdi.range, driver.effect=sco.HDI, density.new=density.new)) #v.weight=v.weight
  #return(list(sir.est=sir.est, obs.vec=obs.vec))
}


