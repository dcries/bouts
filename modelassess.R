#need dgenpois(arma::vec x, arma::vec mu, double lambda, bool logd=true) sourced

pp_assess <- function(mcmcout,Za,nsim,y1real,y2real,weights,burn=1){
  
  ind <- sample(burn:nrow(mcmcout$gamma),nsim)
  n <- nrow(Za)
  nb <- ncol(Za)
  
  #allocate storage
  #y1 posterior predictive checks
  y101 <- rep(0,nsim)
  y110 <- rep(0,nsim)
  y100 <- rep(0,nsim)
  y120 <- rep(0,nsim)
  y102 <- rep(0,nsim)
  y112 <- rep(0,nsim)
  y121 <- rep(0,nsim)
  y122 <- rep(0,nsim)
  y111 <- rep(0,nsim)
  y1meanwpsd <- rep(0,nsim)
  y1wprange <- rep(0,nsim)

  #y2 posterior predictive checks
  kspval <- rep(0,nsim)
  
  #holds compliance rates for raw and weighted data
  pcomply <- rep(0,nsim)
  pcomply2 <- rep(0,nsim)
  
  #holds pi value for each individual during the current iteration
  p <- rep(0,n)
  
  #holds values of distribution of usual MET-minutes for raw and weighted data
  disttable <- matrix(0,nrow=nsim,ncol=39)
  disttablew <- matrix(0,nrow=nsim,ncol=39)
  
  
  #select which values from posterior to use
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    lambda <- mcmcout$lambda[ind]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    for(i in 1:nsim){
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- mvrnorm(n,c(0,0),bcovmat)

      tempmean <- exp(as.numeric((Za%*%gamma[i,]+b[,1])))
      
      x1 <- exp(as.numeric((Za%*%gamma[i,]+b[,1])))
      muy <- Za%*%(betay[i,]) + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Za%*%(betay[i,]) + b[,2]
      muy2 <- Za%*%(betay[i,]) + b[,2]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))

      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
      
      
      x3 <- 30*x1 + p*x2*x1
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum((x3>450/7)*weights)/(sum(weights))
      disttable[i,] <- quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99))
      disttablew[i,] <- wtd.quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99),weights=weights)
      
    }
  

  out <- data.frame(cbind(y1meanwpsd,y1wprange,
                          pcomply,pcomply2,y110,y101,y100,
                          y102,y120,y112,y121,y122,y111))
  return(list(out=out,kspval=kspval,disttable=disttable,disttablew=disttablew))

}