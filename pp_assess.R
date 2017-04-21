#need dgenpois(arma::vec x, arma::vec mu, double lambda, bool logd=true) sourced

pp_assess <- function(mcmcout,Zb,nsim, ymodel,burn=1){
  
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
  
  y2q35 <- rep(0,nsim)
  y2q90 <- rep(0,nsim)
  y2daydiff <- rep(0,nsim)
  #both
  y1y2regcoef <- rep(0,nsim)
  y1y2cor <- rep(0,nsim)
  pcomply <- rep(0,nsim)
  pcomply2 <- rep(0,nsim)
  
  p <- rep(0,n)
  

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
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy+sigma2y[i]/2)
      
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]
      
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
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
    p <- mcmcout$p[ind,]
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
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
    p <- mcmcout$p[ind,]
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
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1] + b[,3]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i]+b[j,2],FALSE,m[i,j])
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]+ b[,3]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]+ b[,3]
      
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
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
    p <- mcmcout$p[ind,]

    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      cov1 = sqrt(sigma2b[i,1]*sigma2b[i,2])*corrb[i]
      bcovmat <- matrix(c(sigma2b[i,1],cov1,cov1,sigma2b[i,2]),ncol=2,byrow=T)
      b <- mvrnorm(n,c(0,0),bcovmat)
      
      tempmean <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #ind <- which(tempmean*(1-lambda[i]+m[i,]*lambda[i]) < 0)
      
      
      x1 <- exp(as.numeric((Zb%*%gamma[i,]+b[,1])))
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1] + b[,2]
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]+ b[,2]
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]+ b[,2]
      
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
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
    p <- mcmcout$p[ind,]
    
    for(i in 1:nsim){
      #corresponds to nci model with correlated re in mu and lambda (normal)
      b <- rnorm(n,0,sqrt(sigma2b[i]))
      
      mux1 <- exp(as.numeric((Zb%*%gamma[i,])))
      x1 <- rgamma(n,eta[i],eta[i]/mux1)
      #print(summary(x1))
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1] + b
      x2 <- exp(muy+sigma2y[i]/2)
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE)
        
      }
      muy1 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1]+ b
      muy2 <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1]+ b
      
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
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
      muy <- Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(x1)*betay[i,ncol(Zb)+1]
      x2 <- exp(muy)
      
      y1 <- matrix(0,ncol=2,nrow=n)
      
      for(j in 1:n){
        #y1[j,] <- rgenpois(2,x1[j],lambda[i]+b[j,2])
        y1[j,] <- rgenpois(2,x1[j],lambda[i])
        p[j] <- 1-dgenpois(0,x1[j],lambda[i],FALSE,Inf)
        
      }
      muy1 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,1])*betay[i,ncol(Zb)+1])
      muy2 <- exp(Zb%*%(betay[i,1:ncol(Zb)]) + sqrt(y1[,2])*betay[i,ncol(Zb)+1])
      
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
      
      y2zeroboth[i] <- sum(rowSums(y2)==0)
      y2zeroeither[i] <- sum(apply(y2,1,function(x){return(!0%in%x)}))
      y2greaterthan[i] <- sum(y1[,1]*30+y2[,1]>450/7)
      y2median[i] <- median(c(y2))
      y2q15[i] <- quantile(c(y2),probs=c(0.20))
      y2q25[i] <- quantile(c(y2),probs=c(0.25))
      y2q35[i] <- quantile(c(y2),probs=c(0.35))
      y2q90[i] <- quantile(c(y2),probs=c(0.9))
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
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
      y2daydiff[i] <- mean(y2[,1]-y2[,2])
      #y1y2regcoef[i] <- coef(lm((y2[,1]-y2[,2])~I(y1[,1]-y1[,2])+weekenddiff))[2]
      y1y2cor[i] <- cor(c(y1[y1>0]),c(y2[y2>0]))
      
      #print(i)
      
      x3 <- 30*x1 + p*x2
      pcomply[i] <- sum(x3>450/7)/n
      pcomply2[i] <- sum(x3>450/5)/n
      
    }
  }
    
  out <- data.frame(cbind(y1zeroboth,y1zeroeither,y1ones,y1twos,y1meanwpsd,y1wprange,y1overallrange,
              y2zeroboth,y2zeroeither,y2greaterthan,y1y2regcoef,y1y2cor,y2median,
              y2q15,y2q25,y2q35,y2q90,y2daydiff,pcomply,pcomply2))
  return(list(out=out,y1=y1,y2=y2,x2=x2,x3=x3,p=p,x1=x1))
  #return(out)
  
}