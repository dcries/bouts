#need dgenpois(arma::vec x, arma::vec mu, double lambda, bool logd=true) sourced

pp_assess <- function(mcmcout,Zb,nsim, ymodel,y1real,y2real,weights,burn=1){
  
  ind <- sample(burn:nrow(mcmcout$gamma),nsim)
  n <- nrow(Zb)
  #na <- ncol(Za)
  nb <- ncol(Zb)
  
  #allocate storage
  #y1 checks
  y1zeroboth <- rep(0,nsim)
  y1zeroeither <- rep(0,nsim)
  y1ones <- rep(0,nsim)
  y1twos <- rep(0,nsim)
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
  y1overallrange <- rep(0,nsim)
  #y2 checks
  y2zeroboth <- rep(0,nsim)
  y2zeroeither <- rep(0,nsim)
  y2greaterthan <- rep(0,nsim)
  y2median <- rep(0,nsim)
  y2q15 <- rep(0,nsim)
  y2q25 <- rep(0,nsim)
  y2q30 <- rep(0,nsim)
  
  y2q35 <- rep(0,nsim)
  y2q90 <- rep(0,nsim)
  y2q95 <- rep(0,nsim)
  y2mean <- rep(0,nsim)
  y2var <- rep(0,nsim)
  #both
  y1y2regcoef <- rep(0,nsim)
  y1y2cor <- rep(0,nsim)
  pcomply <- rep(0,nsim)
  pcomply2 <- rep(0,nsim)
  
  kspval <- rep(0,nsim)
  
  p <- rep(0,n)
  
  disttable <- matrix(0,nrow=nsim,ncol=39)
  disttablew <- matrix(0,nrow=nsim,ncol=39)
  

  if(ymodel==1){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    for(i in 1:nsim){
      #corresponds to nci model with no re for lambda or in y2
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy+sigma2y[i]/2)
      
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]
      
      #check0 <- matrix(rbinom(2*n,1,rep(p,each=2)),ncol=2,byrow=TRUE)
      #y1 <- y1[,1:2]
      check0 <- y1
      check0[check0 > 0] <- 1
      
      # if(ymodel=="gamma"){
      #   y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
      # }
      # else if(ymodel=="lognormal"){
        y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      # }
      # else{
      #   stop("ymodel specification not known")
      # }
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2q95[i] <- quantile(c(y2),probs=c(0.95))
      y2q30[i] <- quantile(c(y2),probs=c(0.3))
      
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
    }
  }
  if(ymodel==2){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    m <- mcmcout$m[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      bcovmat <- matrix(c(sigma2b[i,1],sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i],sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i],sigma2b[i,2]),ncol=2,byrow=T)
      b <- mvrnorm(n,c(0,0),bcovmat)
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      for(j in 1:n){
        while((tempmean[j]*(1-(lambda[i]+b[j,2]))+m[i,j]*(lambda[i]+b[j,2]) < 0) || (lambda[i]+b[j,2] > 1)){
          b[j,] <- mvrnorm(1,c(0,0),bcovmat)
          tempmean[j] <- exp(as.numeric((Zb[j,]%*%gamma[i,]+b[j,1])))
          #print(tempmean[j])
        }
      }
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i]+b[j,2],FALSE,m[i,j])

      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]
      
      check0 <- y1
      check0[check0 > 0] <- 1

      
      # if(ymodel=="gamma"){
      #   y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
      # }
      # else if(ymodel=="lognormal"){
        y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      # }
      # else{
      #   stop("ymodel specification not known")
      # }
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2q95[i] <- quantile(c(y2),probs=c(0.95))
      y2q30[i] <- quantile(c(y2),probs=c(0.3))
      
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
    }
  }
  if(ymodel=="2b"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind,]
    
    #mux2 <- mcmcout$mux2[ind,]
    m <- mcmcout$m[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i,1]
      cov2 = sqrt(sigma2b[i,1]*sigma2b[i,3])*corrb[i,2]
      cov3 = sqrt(sigma2b[i,2]*sigma2b[i,3])*corrb[i,3]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov2,cov1,sigma2b[i,2],cov3,cov2,cov3,sigma2b[i,3]),
                        ncol=3,byrow=T)
      b <- mvrnorm(n,c(0,0,0),bcovmat)
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      for(j in 1:n){
        while((tempmean[j]*(1-(lambda[i]+b[j,2]))+m[i,j]*(lambda[i]+b[j,2]) < 0) || (lambda[i]+b[j,2] > 1)){
          b[j,] <- mvrnorm(1,c(0,0,0),bcovmat)
          tempmean[j] <- exp(as.numeric((Zb[j,]%*%gamma[i,]+b[j,1])))
          #print(tempmean[j])
        }
      }
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,3]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(999,ncol=2,nrow=n)
      
      for(j in 1:n){
        #while(y1[j,1] > 110){
          y1[j,1] <- rgenpois(1,x1[j],lambda[i]+b[j,2])
        #}
        #while(y1[j,2] > 110){
          y1[j,2] <- rgenpois(1,x1[j],lambda[i]+b[j,2])
        #}
        p[j] <- 1-dgenpois(0,x1[j],lambda[i]+b[j,2],FALSE,m[i,j])
        
      }
      #print(summary(lambda[i]+b[,2]))
      
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,3]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,3]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      # if(ymodel=="gamma"){
      #   y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
      # }
      # else if(ymodel=="lognormal"){
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      # }
      # else{
      #   stop("ymodel specification not known")
      # }
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2q95[i] <- quantile(c(y2),probs=c(0.95))
      y2q30[i] <- quantile(c(y2),probs=c(0.3))
      
            #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      #print(i)
      
    }
  }
  if(ymodel=="2c"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]

    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- mvrnorm(n,c(0,0),bcovmat)
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)

      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      

      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)

      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
 
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2q95[i] <- quantile(c(y2),probs=c(0.95))
      y2q30[i] <- quantile(c(y2),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var <- var(y2[y2>0])
      #print(i)

      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="2c2"){
    gamma <- mcmcout$gamma[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      b <- rnorm(n,0,sqrt(sigma2b[i]))
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b)))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      #print(summary(x1))
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      x2 <- rep(0,n)
      x3 <- rep(0,n)
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))      

    }
  }
  if(ymodel=="2c3"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    y1 <- y1real
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      #b <- mvrnorm(n,c(0,0),bcovmat)
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      # x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      # #print(summary(x1))
      # muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      # x2 <- exp(muy+sigma2y[i]/2)
      # y1 <- matrix(0,ncol=2,nrow=n)
      # 
      # for(j in 1:n){
      #   y1[j,] <- rgenpois(2,x1[j],lambda[i])
      #   p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
      #   
      # }
      b <- rnorm(n,0,sqrt(sigma2b[i]))
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b
      
      #check0 <- y1
      #check0[check0 > 0] <- 1
      muy1 <- muy1[exp(muy1)>0]
      muy2 <- muy2[exp(muy2)>0]
      nm <- length(muy1)+length(muy2)
      #print(nm)
      #y2 <- rlnorm(nm,c(muy1,muy2),rep(sqrt(sigma2y[i]),nm))
      y2 <- rgamma(nm,rep(delta[i],nm),rep(delta[i],nm)/c(exp(muy1),exp(muy2)))
      # y1zeroboth[i] <- sum(rowSums(y1)==0)
      # y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      # y1ones[i] <- sum(y1==1)
      # y1twos[i] <- sum(y1==2)
      # y1meanwpsd[i] <- mean(apply(y1,1,sd))
      # y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      # y1overallrange[i] <- max(y1)-min(y1)
      # y110[i] <- sum(apply(y1,1,function(x){return(x[1]>=1 & x[2]==0)}))
      # y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=1)}))
      
      # y2zeroboth[i] <- sum(rowSums(y2)==0)
      # y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      #y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
      
      #print(i)
      x1 <- rep(0,n)
      x2 <- rep(0,n)
      x3 <- 30*x1 + x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
    }
  }
  if(ymodel=="2d"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind]

    #mux2 <- mcmcout$mux2[ind,]

    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      b <- rnorm(n,0,sqrt(sigma2b[i]))
      
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      # if(ymodel=="gamma"){
      #   y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
      # }
      # else if(ymodel=="lognormal"){
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      # }
      # else{
      #   stop("ymodel specification not known")
      # }
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
    }
  }
  
  if(ymodel==3){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    

    for(i in 1:nsim){
      #corresponds to nci model with no re for lambda or in y2
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy)
      
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1])
      muy2 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1])
      
      #check0 <- matrix(rbinom(2*n,1,rep(p,each=2)),ncol=2,byrow=TRUE)
      #y1 <- y1[,1:2]
      check0 <- y1
      check0[check0 > 0] <- 1
      
      # if(ymodel=="gamma"){
      #   y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),rep(delta[i],2*n)/rep(x2,each=2)),ncol=2,byrow=TRUE)
      # }
      # else if(ymodel=="lognormal"){
      y2 <- check0*matrix(rgamma(2*n,rep(delta[i],2*n),delta[i]/c(muy1,muy2)),ncol=2,byrow=FALSE)
      # }
      # else{
      #   stop("ymodel specification not known")
      # }
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
      
    }
  }
  if(ymodel==4){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]

    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    m <- mcmcout$m[ind,]

    for(i in 1:nsim){
      #corresponds to nci model with  re for lambda 
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy+sigma2y[i]/2)
      b <- rnorm(n,0,sqrt(sigma2b[i]))

      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i]+b[j],FALSE,m[i,j])
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]

      #check0 <- matrix(rbinom(2*n,1,rep(p,each=2)),ncol=2,byrow=TRUE)
      #y1 <- y1[,1:2]
      check0 <- y1
      check0[check0 > 0] <- 1

      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)


      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2q95[i] <- quantile(c(y2),probs=c(0.95))
      y2q30[i] <- quantile(c(y2),probs=c(0.3))
      
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
    }
  }
  if(ymodel=="5"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- matrix(0,nrow=n,ncol=2)
      b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      b[which(b[,1]>6),1] <- 6
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
 
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="5b"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- matrix(0,nrow=n,ncol=2)
      b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- (muy+sigma2y[i]/2)^2
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- log((Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2])^2)
      muy2 <- log((Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2])^2)
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0

      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="6"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- matrix(0,nrow=n,ncol=2)
      b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      #b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- exp(muy)*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2])
      muy2 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2])
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- check0*matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="7"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- matrix(0,nrow=n,ncol=2)
      b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      b[which(b[,1]>6),1] <- 6
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #print(x1[j]);print(b[j,1])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- (Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2])^(2)
      muy2 <- (Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2])^(2)
      check0 <- y1
      check0[check0 > 0] <- 1
      
      y2 <- matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum((x3>450/7)*weights)/sum(weights)
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="7b"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- matrix(0,nrow=n,ncol=2)
      b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- (Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2])^(2)
      muy2 <- (Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2])^(2)
      check0 <- y1
      check0[check0 > 0] <- 1
      
      y2 <- matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="7c"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- mvrnorm(n,c(0,0),bcovmat)
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- (Zb%*%(betay[i,1:ncol(Zb)]) + b[,2])^(2)
      muy2 <- (Zb%*%(betay[i,1:ncol(Zb)]) + b[,2])^(2)
      check0 <- y1
      check0[check0 > 0] <- 1
      
      y2 <- matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum((x3>450/7)*weights)/sum(weights)
      
      disttable[i,] <- quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99))
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="8"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      #b <- matrix(0,nrow=n,ncol=2)
      #b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      #[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1]
      
      
      #x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      #muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- (Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,1])*betay[i,ncol(Zb)+1])^(2)
      muy2 <- (Zb%*%(betay[i,1:ncol(Zb)]) + log(y1[,2])*betay[i,ncol(Zb)+1])^(2)
      check0 <- y1
      check0[check0 > 0] <- 1
      
      y2 <- matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="8b"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    phi <- mcmcout$phi[ind]
    
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      #cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      #bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      #b <- matrix(0,nrow=n,ncol=2)
      #b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      #[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      
      #tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1]
      
      
      #x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      #muy <- Zb%*%(betay[i,1:ncol(Zb)]) + log(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- ((muy)^(2))*gamma(1+1/phi[i])
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- (Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1])^(2)
      muy2 <- (Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1])^(2)
      check0 <- y1
      check0[check0 > 0] <- 1
      
      y2 <- matrix(rweibull(2*n,rep(phi[i],2*n),c(muy1,muy2)),ncol=2,byrow=FALSE)
      y2[check0==0] <- 0
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
    }
  }
  if(ymodel=="9"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      #b <- matrix(0,nrow=n,ncol=2)
      b <- mvrnorm(n,c(0,0),bcovmat)
      #b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      #b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      #b[which(b[,1]>6),1] <- 6
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,]) + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #print(x1[j]);print(lambda[i])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,]) + b[,2]
      muy2 <- Zb%*%(betay[i,]) + b[,2]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2*x1
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum((x3>450/7)*weights)/(sum(weights))
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
      disttable[i,] <- quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99))
      disttablew[i,] <- wtd.quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99),weights=weights)
      
    }
  }
  if(ymodel=="9b"){
    gamma <- mcmcout$gamma[ind,]
    betay <- mcmcout$betay[ind,]
    alpha <- mcmcout$alpha[ind,]
    sigma2y <- mcmcout$sigma2y[ind]
    sigma2x <- mcmcout$sigma2x[ind]
    
    eta <- mcmcout$eta[ind]
    theta <- mcmcout$theta[ind]
    delta <- mcmcout$delta[ind]
    lambda <- mcmcout$lambda[ind]
    latentx1 <- mcmcout$latentx1[ind,]
    latentx2 <- mcmcout$latentx2[ind,]
    #mux1 <- mcmcout$mux1[ind,]
    #muy <- mcmcout$muy[ind,]
    sigma2b <- mcmcout$sigma2b[ind,]
    corrb <- mcmcout$corrb[ind]
    
    #mux2 <- mcmcout$mux2[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      #b <- matrix(0,nrow=n,ncol=2)
      b <- mvrnorm(n,c(0,0),bcovmat)
      #b[,1] <- rnorm(n,0,sqrt(sigma2b[i,1]))
      #b[,2] <- rnorm(n,0,sqrt(sigma2b[i,2]))
      #b[which(b[,1]>6),1] <- 6
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,]) + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #print(x1[j]);print(lambda[i])
        y1[j,] <- rpois(2,x1[j])
        p[j] <- 1-dpois(0,x1[j],TRUE)
        
      }
      muy1 <- Zb%*%(betay[i,]) + b[,2]
      muy2 <- Zb%*%(betay[i,]) + b[,2]
      
      check0 <- y1
      check0[check0 > 0] <- 1
      
      
      y2 <- check0*matrix(rlnorm(2*n,c(muy1,muy2),rep(sqrt(sigma2y[i]),2*n)),ncol=2,byrow=FALSE)
      
      y1zeroboth[i] <- sum(rowSums(y1)==0)
      y1zeroeither[i] <- sum(apply(y1,1,function(x){return(!0%in%x)}))
      y1ones[i] <- sum(y1==1)
      y1twos[i] <- sum(y1==2)
      y1meanwpsd[i] <- mean(apply(y1,1,sd))
      y1wprange[i] <- max(abs(apply(y1,1,function(x){return(x[2]-x[1])})))
      y1overallrange[i] <- max(y1)-min(y1)
      
      y100[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==0)}))
      y110[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==0)}))
      y101[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]==1)}))
      y120[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==0)}))
      y102[i] <- sum(apply(y1,1,function(x){return(x[1]==0 & x[2]>=2)}))
      y112[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]>=2)}))
      y121[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]==1)}))
      y111[i] <- sum(apply(y1,1,function(x){return(x[1]==1 & x[2]==1)}))
      y122[i] <- sum(apply(y1,1,function(x){return(x[1]>=2 & x[2]>=2)}))
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2[y2>0]))
      y2q15[i] <- quantile(c(y2[y2>0]),probs=c(0.20))
      y2q25[i] <- quantile(c(y2[y2>0]),probs=c(0.25))
      y2q35[i] <- quantile(c(y2[y2>0]),probs=c(0.35))
      y2q90[i] <- quantile(c(y2[y2>0]),probs=c(0.9))
      y2q95[i] <- quantile(c(y2[y2>0]),probs=c(0.95))
      y2q30[i] <- quantile(c(y2[y2>0]),probs=c(0.3))
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      y2mean[i] <- mean(y2[y2>0])
      y2var[i] <- var(y2[y2>0])
      #print(i)
      
      x3 <- 30*x1 + p*x2*x1
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum((x3>450/7)*weights)/(sum(weights))
      kspval[i] <- ks.test(y2real[y2real>0],y2[y2>0])$p.value
      disttable[i,] <- quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99))
      disttablew[i,] <- wtd.quantile(x3,probs=c(0.01,seq(from=0.05,to=0.95,by=0.025),0.99),weights=weights)
      
    }
  }
  out <- data.frame(cbind(y1zeroboth,y1zeroeither,y1ones,y1twos,y1meanwpsd,y1wprange,y1overallrange,
              y2zeroboth,y2zeroeither,y2greaterthan,y1y2regcoef,y1y2cor,y2median,
              y2q15,y2q25,y2q35,y2q90,y2q95,y2q30,pcomply,pcomply2,y110,y101,y100,
              y102,y120,y112,y121,y122,y111,y2mean,y2var))
  return(list(out=out,y1=y1,y2=y2,x2=x2,x3=x3,p=p,x1=x1,kspval=kspval,disttable=disttable,disttablew=disttablew))
  #return(out)
  
}