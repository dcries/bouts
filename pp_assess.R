#need dgenpois(arma::vec x, arma::vec mu, double lambda, bool logd=true) sourced

pp_assess <- function(mcmcout,Zb,weekenddiff,nsim, ymodel){
  
  ind <- sample(1:nrow(mcmcout$latentx1),nsim)
  
  gamma <- mcmcout$gamma[ind,]
  betax <- mcmcout$betax[ind,]
  alpha <- mcmcout$alpha[ind,]
  sigma2y <- mcmcout$sigma2y[ind]
  sigma2x <- mcmcout$sigma2x[ind]
  
  eta <- mcmcout$eta[ind]
  theta <- mcmcout$theta[ind]
  delta <- mcmcout$delta[ind]
  lambda <- mcmcout$lambda[ind]
  latentx1 <- mcmcout$latentx1[ind,]
  latentx2 <- mcmcout$latentx2[ind,]
  mux1 <- mcmcout$mux1[ind,]
  #mux2 <- mcmcout$mux2[ind,]
  p <- mcmcout$p[ind,]
  
  n <- ncol(latentx1)
  #na <- ncol(Za)
  nb <- ncol(Zb)
  
  #allocate storage
  #y1 checks
  y1zeroboth <- rep(0,nsim)
  y1zeroeither <- rep(0,nsim)
  y1meanwpsd <- rep(0,nsim)
  y1wprange <- rep(0,nsim)
  y1overallrange <- rep(0,nsim)
  #y2 checks
  y2zeroboth <- rep(0,nsim)
  y2zeroeither <- rep(0,nsim)
  y2greaterthan <- rep(0,nsim)
  y2median <- rep(0,nsim)
  y2q35 <- rep(0,nsim)
  y2q90 <- rep(0,nsim)
  #both
  y1y2regcoef <- rep(0,nsim)
  pcomply <- rep(0,nsim)
  pcomply2 <- rep(0,nsim)
  
  p2 <- rep(0,nsim)
  
  for(i in 1:nsim){
    #p <- pnorm(Zb%*%(alpha[i,]),0,1)
    #pmat <- matrix(rep(p,each=2),ncol=2,byrow=TRUE)
    x1 <- rgamma(n,eta[i],eta[i]/mux1[i,])
    p <- 1-dgenpois(rep(0,n),x1,lambda[i],FALSE)
    
    mux2 <- exp(Zb%*%betax[i,-(nb+1)] + betax[i,nb+1]*x1)
    x2 <- rgamma(n,theta[i],theta[i]/mux2)
    #x2 <- p*rlnorm(n,log(mux2),sqrt(sigma2x))
    
    y1 <- matrix(0,ncol=2,nrow=n)
    
    for(j in 1:n){
      y1[j,] <- rgenpois(2,x1[j],lambda[i])
    }
    
    #check0 <- matrix(rbinom(2*n,1,rep(p,each=2)),ncol=2,byrow=TRUE)
    check0 <- y1
    check0[check0 > 0] <- 1
    
    if(ymodel=="gamma"){
      y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
    }
    else if(ymodel=="lognormal"){
      y2 <- check0*matrix(rlnorm(2*n,rep(log(x2),each=2),rep(sqrt(sigma2y),2*n)),ncol=2,byrow=TRUE)
    }
    else{
      stop("ymodel specification not known")
    }
    
    y1zeroboth[i] <- sum(rowSums(y1)==0)
    y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
    y1meanwpsd[i] <- mean(apply(y1,1,sd))
    y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
    y1overallrange[i] <- max(y1)-min(y1)
    
    y2zeroboth[i] <- sum(rowSums(y2)==0)
    y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
    y2greaterthan[i] <- sum(y2>450/7)
    y2median[i] <- median(c(y2))
    y2q35[i] <- quantile(c(y2),probs=c(0.35))
    y2q90[i] <- quantile(c(y2),probs=c(0.9))
    
    y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
    #print(i)
    
    x3 <- 30*x1 + p*x2
    pcomply[i] <- sum(x3>450/7)/n
    pcomply2[i] <- sum(x3>450/5)/n
    
    #-----------------------------------
    # p1 <- matrix(1,nrow=n,ncol=7)
    # for(j in 1:7){
    #   x1 <- rgamma(n,eta[i],eta[i]/mux1[i,])
    #   mux2 <- exp(Zb%*%betax[i,-(nb+1)] + betax[i,nb+1]*x1)
    #   
    #   #p <- pnorm(Zb%*%(alpha[i,]),0,1)
    #   p <- 1-dgenpois(rep(0,n),latentx1[i,],lambda[i],FALSE)
    #   x2 <- rgamma(n,theta[i],theta[i]/mux2)
    #   #x2 <- rlnorm(n,log(mux2),sqrt(sigma2x))
    #   x3 <- 30*x1 + p*x2
    #   p1[x3 < 90,j] <- 0
    # }
    # p2[i] <- sum(rowSums(p1)>=5)/nrow(p1)
    # 
  }
  out <- data.frame(cbind(y1zeroboth,y1zeroeither,y1meanwpsd,y1wprange,y1overallrange,
              y2zeroboth,y2zeroeither,y2greaterthan,y1y2regcoef,y2median,
              y2q35,y2q90,pcomply,pcomply2))
  return(list(out=out,y1=y1,y2=y2,x2=x2,mux2=mux2,x3=x3,p=p,x1=x1))
  #return(out)
  
}