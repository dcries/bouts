

pp_assess <- function(mcmcout,Zb,weekenddiff,nsim){
  
  ind <- sample(1:nrow(mcmcout$latentx1),nsim)
  
  gamma <- mcmcout$gamma[ind,]
  betax <- mcmcout$betax[ind,]
  alpha <- mcmcout$alpha[ind,]
  sigma2y <- mcmcout$sigma2y[ind]
  eta <- mcmcout$eta[ind]
  theta <- mcmcout$theta[ind]
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
  #both
  y1y2regcoef <- rep(0,nsim)
  pcomply <- rep(0,nsim)
  
  for(i in 1:nsim){
    x1 <- rgamma(n,eta[i],eta[i]/mux1[i,])
    mux2 <- exp(Zb%*%betax[i,-(nb+1)] + betax[i,nb+1]*x1)
    x2 <- rgamma(n,theta[i],theta[i]/mux2)
    
    y1 <- matrix(0,ncol=2,nrow=n)
    
    for(j in 1:n){
      y1[j,] <- rgenpois(2,x1[j],lambda[i])
    }
    
    #p <- pnorm(alpha[i,1]+alpha[i,2]*x1+alpha[i,3]*x2,0,1)
    #check0 <- matrix(rbinom(2*n,1,rep(p,each=2)),ncol=2,byrow=TRUE)
    check0 <- y1
    check0[check0 > 0] <- 1
    
    y2 <- check0*matrix(rlnorm(2*n,rep(log(x2),each=2),rep(sqrt(sigma2y),2*n)),ncol=2,byrow=TRUE)
    
    y1zeroboth[i] <- sum(rowSums(y1)==0)
    y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
    y1meanwpsd[i] <- mean(apply(y1,1,sd))
    y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
    y1overallrange[i] <- max(y1)-min(y1)
    
    y2zeroboth[i] <- sum(rowSums(y2)==0)
    y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
    y2greaterthan[i] <- sum(y2>450/7)
    
    y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
    #print(i)
    
    x3 <- 30*x1 + x2
    pcomply[i] <- sum(x3>450/7)/n
  }
  out <- data.frame(cbind(y1zeroboth,y1zeroeither,y1meanwpsd,y1wprange,y1overallrange,
              y2zeroboth,y2zeroeither,y2greaterthan,y1y2regcoef,pcomply))
  return(list(out=out,y1=y1,y2=y2,x2=x2,mux2=mux2,x3=x3))
  return(out)
  
}