library(MASS)
library(MCMCpack)
library(mvtnorm)

pams <- read.csv("C:\\Users\\dcries\\github\\pams\\FinalPAMSDataSetNew.csv")
pams <- pams[-which(is.na(pams$ModPAR)|is.na(pams$Age)|is.na(pams$BMI)|is.na(pams$Smoker)),]


log_qu <- function(y,b,p,mu,sigma2,Sigmab){
  ll1 <- dmvnorm(b,c(0,0),Sigmab,log=TRUE)
  if(y==0){
    ll2 <- log(1-p)
  }
  else{
    ll2 <- log(p) + dlnorm(y,mu,sqrt(sigma2),log=TRUE)
  }
  return(ll1+ll2)
}

log_qx <- function(x,w,y,p,mu,gamma,sigma2,sigma2w,sigma2y){ 
  if(x==0){
    ll1 <- log(1-p)
  }
  else{
    ll1 <- log(p) + dlnorm(x,mu,sqrt(sigma2),log=TRUE)
  } 
  ll3 <- dnorm(w,x,sqrt(sigma2w),log=TRUE) + dnorm(y,x*gamma,sqrt(sigma2y),log=TRUE)
  
  return(ll1+ll3)
}  

log_gx <- function(x,p0,scl){
  if(x==0){
    ll <- log(1-p0)
  }
  else{
    ll <- log(p0) + 2*dcauchy(x,scale=scl,log=TRUE)
  }
  return(ll)
}

mcmc_ln <- function(data,init,prior,nreps,burn=1000){
  #!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!
  #non adaptive right now
  p0 <- 0.15
  
  #data
  Xa    <- data$Xa
  Xb <- data$Xb
  y <- data$y
  w <- data$w
  n <- nrow(Xa)
  z <- rep(1,n)
  z[y==0] <- 0
  ind <- y>0
  yp <- y[ind]
  Xbp <- Xb[ind,]
  nr <- ncol(y)
  ybar <- rowMeans(y)
  scl <- quantile(rowMeans(w),probs=0.75)
  
  #initial values
  currentbeta <- init$currentbeta
  currentalpha <- init$currentalpha
  currentsigma2 <- init$currentsigma2
  currentSigmab <- init$currentSigmab
  currentb <- init$currentb
  currentsigma2w <- init$currentsigma2w
  currentsigma2y <- init$currentsigma2y
  currentgamma <- init$currentgamma
  currentx <- init$currentx
  currentu <- rep(0,n)
  
  #priors
  mu0a <- prior$mu0a
  mu0b <- prior$mu0b
  mu0g <- prior$mu0g
  V0a <- prior$V0a
  V0b <- prior$V0b
  V0g <- prior$V0g
  a0 <- prior$a0
  b0 <- prior$b0
  a0w <- prior$a0w
  b0w <- prior$b0w
  a0y <- prior$a0y
  b0y <- prior$b0y
  nu0 <- prior$nu0
  D0 <- prior$D0
  
  #allocate storae
  beta <- matrix(0,ncol=ncol(Xb),nrow=nreps)
  alpha <- matrix(0,ncol=ncol(Xa),nrow=nreps)
  sigma2 <- rep(0,nreps)
  #u <- matrix(0,)
  sigma2b1 <- rep(0,nreps)
  sigma2b2 <- rep(0,nreps)
  corrb <- rep(0,nreps)
  b1 <- matrix(0,ncol=n,nrow=nreps)
  b2 <- matrix(0,ncol=n,nrow=nreps)
  gamma <- rep(0,nreps)
  sigma2w <- rep(0,nreps)
  sigma2y <- rep(0,nreps)
  latentx <- matrix(0,ncol=n,nrow=nreps)
  
  
  tune <- array(0.01*diag(2),dim=c(2,2,n))
  
  for(i in 1:nreps){
    for(j in 1:n){
      rn <- rnorm(1,Xa[j,]%*%currentalpha+currentb[j,1],1)
      currentu[j] <- ifelse(currentx[j]==0,-abs(rn),abs(rn))
    }
    
    #sample alphas
    Va <- solve(solve(V0a)+t(Xa)%*%Xa)
    mua <- Va%*%(solve(V0a)%*%mu0a + t(Xa)%*%(currentu-currentb[,1]))
    currentalpha <- mvrnorm(1,mua,Va)
    currentp <- pnorm(Xa%*%currentalpha+currentb[,1])
    
    #sample betas
    Vb <- solve(solve(V0b)+t(Xbp)%*%Xbp/currentsigma2)
    mub <- Vb%*%(solve(V0b)%*%mu0b + t(Xbp)%*%(log(yp)-currentb[ind,2])/currentsigma2)
    currentbeta <- mvrnorm(1,mub,Vb)
    currentmu <- exp(Xb%*%currentbeta+currentb[,2])
    
    #sample sigma2
    currentsigma2 <- rinvgamma(1,n/2+a0,(b0+0.5*sum((log(yp)-Xbp%*%currentbeta-currentb[ind,2])^2)))
    
    
    #sample covaraicne matrix for re
    currentSigmab <- riwish(nu0+n,D0+t(currentb)%*%currentb)
    
    for(j in 1:n){
      propb <- mvrnorm(1,currentb[j,],2.88*tune[,,j])
      lacceptprobu <- log_qu(currentx[j],propb,currentp[j],currentmu[j],currentsigma2,currentSigmab) - 
        log_qu(currentx[j],currentb[j,],currentp[j],currentmu[j],currentsigma2,currentSigmab)
      
      if(lacceptprobu > log(runif(1))){
        currentb[j,] <- propb
      }
      if((i>20)&(i<burn)&(i%%20==0)){
        tune[,,j] <- cov(cbind(b1[1:(i-1),j],b2[1:(i-1),j]))
      }
      
    }
    
    #sample sigma2y sigma2w
    wsum <- 0
    ysum <- 0
    for(j in 1:nr){
      wsum <- wsum + sum((w[,j]-currentx)^2)
      ysum <- ysum + sum((y[,j]-currentgamma*currentx)^2)
    }
    currentsigma2w <- rinvgamma(1,n/2+a0w,b0w+0.5*wsum)
    currentsigma2y <- rinvgamma(1,n/2+a0y,b0y+0.5*ysum)
    
    #sample gamma
    Vg <- solve(t(currentx)%*%currentx/currentsigma2y + 1/V0g)
    mg <- Vg*(mu0g/V0g+t(currentx)%*%ybar/currentsigma2y)
    currentgamma <- rnorm(1,mg,Vg)
    
    for(j in 1:n){
      r <- runif(1)
      propx <- ifelse(r<(1-currentp[j]),0,abs(rcauchy(1,scale=scl)))
      lacceptprobx <- log_qx(propx,w[j,],y[j,],currentp[j],currentmu[j],currentgamma,currentsigma2,currentsigma2w,currentsigma2y) + log_gx(currentx[j],p0,scl) -
                      log_qx(currentx[j],w[j,],y[j,],currentp[j],currentmu[j],currentgamma,currentsigma2,currentsigma2w,currentsigma2y,scl,p0) - log_gx(propx,p0,scl)
      
      if(lacceptprobx > log(runif(1))){
        currentx[j] <- propx
      }
    }
    
    
    
    beta[i,] <- currentbeta
    alpha[i,] <- currentalpha
    sigma2[i] <- currentsigma2
    #u <- matrix(0,)
    sigma2b1[i] <- currentSigmab[1,1]
    sigma2b2[i] <- currentSigmab[2,2]
    corrb[i] <- currentSigmab[1,2]/sqrt(currentSigmab[1,1]*currentSigmab[2,2])
    b1[i,] <- t(currentb[,1])
    b2[i,] <- t(currentb[,2])
    sigma2w[i] <- currentsigma2w
    sigma2y[i] <- currentsigma2y
    gamma[i] <- currentgamma
    latentx[i,] <- currentx
    
    if(i%%100==0){
      print(i)
    }
  }
  params <- data.frame(cbind(alpha,beta,sigma2,sigma2b1,sigma2b2,corrb))
  names(params) <- c(paste0("alpha_",0:(ncol(Xa)-1)),
                     paste0("beta_",0:(ncol(Xb)-1)),
                     "sigma2","sigma2b1","sigma2b2","corrb")
  
  out <- list(params=params,b1=b1,b2=b2)
  return(out)
} 

#pams <- pams[1:500,]
data <- list(Xa=with(pams,model.matrix(~Age+BMI+Gender+Smoker)),
             Xb=with(pams,model.matrix(~Age+BMI+Gender+Smoker)),
             y=pams$ModPAR)
init <- list(currentbeta=rep(0,ncol(data$Xb)),
             currentalpha=rep(0,ncol(data$Xa)),
             currentsigma2=1,
             currentSigmab=diag(2),
             currentb=matrix(0,ncol=2,nrow=nrow(pams)))
prior <- list(mu0a=rep(0,ncol(data$Xa)),
              mu0b=rep(0,ncol(data$Xb)),
              V0a=10*diag(ncol(data$Xa)),
              V0b=10*diag(ncol(data$Xb)),
              a0=1,
              b0=1,
              nu0=3,
              D0=diag(2))

out <- mcmc_ln(data,init,prior,10000)
