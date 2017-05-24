

usual_ee <- function(mcmcout,Zorig,nsim){
  ind <- sample(burn:nrow(mcmcout$gamma),nsim)
  
  gamma <- mcmcout$gamma[ind,]
  betay <- mcmcout$betay[ind,]
  phi <- mcmcout$phi[ind]
  lambda <- mcmcout$lambda[ind]
  sigma2b <- mcmcout$sigma2b[ind,]
  corrb <- mcmcout$corrb[ind]
  n <- nrow(Zorig)
  
  pcomply <- matrix(0,ncol=18,nrow=nsim) #18 different groups
  disttable <- array(0,c(nsim,39,18))
  
  for(i in 1:nsim){
    #simulate covariate matricies for different groups
    Zmale <- Zorig; Zmale[,3] <- rep(1,n)
    Zfemale <- Zorig; Zfemale[,3] <- rep(0,n)
    Z20 <- Zorig; Z20[,2] <- sample(20:29,n,replace=TRUE)
    Z30 <- Zorig; Z30[,2] <- sample(30:39,n,replace=TRUE)
    Z40 <- Zorig; Z40[,2] <- sample(40:49,n,replace=TRUE)
    Z50 <- Zorig; Z50[,2] <- sample(50:59,n,replace=TRUE)
    Z60 <- Zorig; Z60[,2] <- sample(60:70,n,replace=TRUE)
    Zbmi25 <- Zorig; Zbmi25[,4] <- sample(Zorig[Zorig[,4]<25,4],n,replace=TRUE)
    Zbmi30 <- Zorig; Zbmi30[,4] <- sample(Zorig[(Zorig[,4]>25) && (Zorig[Zorig[,4]<30]),4],n,replace=TRUE)
    Zbmi35 <- Zorig; Zbmi35[,4] <- sample(Zorig[Zorig[,4]>30 && Zorig[Zorig[,4]<30],4],n,replace=TRUE)
    Zbmi35p <- Zorig; Zbmi35p[,4] <- sample(Zorig[Zorig[,4]>35,4],n,replace=TRUE)
    Zwhite <- Zorig; Zwhite[,7:8] <- 0
    Zblack <- Zorig; Zblack[,7] <- 1; Zblack[,8] <- 0
    Zhispanic <- Zorig; Zhispanic[,7] <- 0; Zhispanic[,8] <- 1
    Zcollege <- Zorig; Zcollege[,6] <- 1
    Znocollege <- Zorig; Znocollege[,6] <- 0
    Zsmoke <- Zorig; Zsmoke[,5] <- 1
    Znosmoke <- Zorig; Znosmoke[,5] <- 0
    
    Z <- list(Zmale,Zfemale,Z20,Z30,Z40,Z50,Z60,Zbmi25,Zbmi30,Zbmi35,Zbmi35p,Zwhite,Zblack,Zhispanic,Zcollege,
              Znocollege,Zsmoke,Znosmoke)
    m <- length(Z)
    #corresponds to nci model with correlated re in mu and lambda (normal)
    cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
    bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
    b <- mvrnorm(n,c(0,0),bcovmat)
  
    for(j in 1:m){
      x1 <- exp(as.numeric((Z[[j]]%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Z[[j]]%*%(betay[i,1:ncol(Zb)]) + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      
      for(k in 1:n){
        p[k] <- 1-dgenpois(0,x1[k],lambda[k],FALSE)
      }
      
      x3 <- 30*x1 + p*x2
      pcomply[i,j] <- sum(x3>450/7)/n
      disttable[i,,j] <- quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99))
      
    }
  }
  return(list(pcomply=pcomply,disttable=disttable))
}