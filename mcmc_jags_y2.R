library(rjags)

modeltp <- "model{
  
  # For the ones trick
  C <- 10000
  
  # for every observation
  for(i in 1:N){
    
    logit(w[i]) <- zeta[i]
    zeta[i] <- alpha0 + alpha1*x1[i] + alpha2*x2[i]
    
    eta[i] <- beta0 + beta1*age[i] + beta2*bmi[i] + beta3*gender[i] + beta4*smoke[i] + beta5*x1[i]
    muy[i] <- theta0 + theta1*x1[i] + theta2*x2[i]

    mux1[i] <- gamma0 + gamma1*age[i] + gamma2*bmi[i] + gamma3*gender[i] + gamma4*smoke[i]
    rate[i] <- shape/exp(mux1[i])

    x2[i] ~ dlnorm(eta[i],tau2x)
    x1[i] ~ dgamma(shape,rate[i])

  }

  for(i in 1:2*N){
      logNormal[i] <- log(dlnorm(y2[i],muy[ind[i]],tau2y))

      logLik[i] <- (1 - z[i]) * log(1 - w[ind[i]]) + z[i] * ( log(w[ind[i]]) + logNormal[i] )
      Lik[i] <- exp(logLik[i])
    
      # Use the ones trick
      p[i] <- Lik[i] / C
      ones[i] ~ dbern(p[i])

      y1[i] ~ dpois(x1[ind[i]])

  }


  
  # PRIORS
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
  beta2 ~ dnorm(0, 0.0001)
  beta3 ~ dnorm(0, 0.0001)
  beta4 ~ dnorm(0, 0.0001)
  beta5 ~ dnorm(0, 0.0001)


  gamma0 ~ dnorm(0, 0.0001)
  gamma1 ~ dnorm(0, 0.0001)
  gamma2 ~ dnorm(0, 0.0001)
  gamma3 ~ dnorm(0, 0.0001)
  gamma4 ~ dnorm(0, 0.0001)

  tau2x ~ dgamma(1, 1)
  tau2y ~ dgamma(1, 1)

  theta0 ~ dnorm(0, 0.0001)
  theta1 ~ dnorm(0, 0.0001)
  theta2 ~ dnorm(0, 0.0001)

  alpha0 ~ dnorm(0, 0.0001)
  alpha1 ~ dnorm(0, 0.0001)
  alpha2 ~ dnorm(0, 0.0001)

  shape ~ dgamma(1,1)

}"


z1 <- rep(1,nrow(data$y2))
z2 <- rep(1,nrow(data$y2))
z1[data$y2[,1]==0] <- 0
z2[data$y2[,2]==0] <- 0
z <- c(z1,z2)

dattp <- list(y1=c(data$y1[,1],data$y1[,2]),y2=c(data$y2[,1],data$y2[,2]),z=z,age=Za[,2],gender=Za[,3],bmi=Za[,4],smoke=Za[,5],
              N=nrow(data$y2),ones=rep(1,2*nrow(data$y2)),ind=rep(1:nrow(data$y2),2))

mtp = jags.model(textConnection(modeltp), dattp,n.adapt=1000,n.chains=3)
rtp = coda.samples(mtp, c("shape","beta0","beta1","beta2","beta3","beta4","tau2y","gamma0","gamma1","gamma2","gamma3","gamma4","tau2x","theta0","theta1"), n.iter=1000)


#-----------------------------------
#no measurement error part
#------------------------------------


modeltp2 <- "model{
  
# For the ones trick
C <- 10000

# for every observation
for(i in 1:N){

# define the logistic regression model, where w is the probability of occurance.
# use the logistic transformation exp(z)/(1 + exp(z)), where z is a linear function
logit(w[i]) <- zeta[i]
zeta[i] <- gamma0 + gamma1*age[i] + gamma2*bmi[i] + gamma3*gender[i] + gamma4*smoke[i]

eta[i] <- beta0 + beta1*age[i] + beta2*bmi[i] + beta3*gender[i] + beta4*smoke[i]


logNormal[i] <- log(dlnorm(y[i],eta[i],tau2x))


logLik[i] <- (1 - z[i]) * log(1 - w[i]) + z[i] * ( log(w[i]) + logNormal[i] )
Lik[i] <- exp(logLik[i])

# Use the ones trick
p[i] <- Lik[i] / C
ones[i] ~ dbern(p[i])

}



# PRIORS
beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)
beta3 ~ dnorm(0, 0.0001)
beta4 ~ dnorm(0, 0.0001)

gamma0 ~ dnorm(0, 0.0001)
gamma1 ~ dnorm(0, 0.0001)
gamma2 ~ dnorm(0, 0.0001)
gamma3 ~ dnorm(0, 0.0001)
gamma4 ~ dnorm(0, 0.0001)

tau2x ~ dgamma(1, 1)


}"

dattp2 <- list(y=c(data$y2[,1],data$y2[,2]),z=z,age=rep(Za[,2],2),gender=rep(Za[,3],2),bmi=rep(Za[,4],2),smoke=rep(Za[,5],2),
              N=2*nrow(data$y2),ones=rep(1,2*nrow(data$y2)))

mtp2 = jags.model(textConnection(modeltp2), dattp2,n.adapt=1000,n.chains=3)
rtp2 = coda.samples(mtp2, c("beta0","beta1","beta2","beta3","beta4","gamma0","gamma1","gamma2","gamma3","gamma4","tau2x"), n.iter=1000)

#-------------------------------

modeltp <- "model{
  
# For the ones trick
C <- 10000

# for every observation
for(i in 1:N){

logit(w[i]) <- zeta[i]
zeta[i] <- alpha0  + alpha2*x2[i]

eta[i] <- beta0 + beta1*age[i] + beta2*bmi[i] + beta3*gender[i] + beta4*smoke[i] 
muy[i] <- theta0  + theta2*x2[i]

#mux1[i] <- gamma0 + gamma1*age[i] + gamma2*bmi[i] + gamma3*gender[i] + gamma4*smoke[i]
#rate[i] <- shape/exp(mux1[i])

rate2[i] <- shape2/exp(eta[i])

x2[i] ~ dgamma(shape2,rate2[i])
#x1[i] ~ dgamma(shape,rate[i])

}

for(i in 1:2*N){
logNormal[i] <- log(dlnorm(y2[i],muy[ind[i]],tau2y))

logLik[i] <- (1 - z[i]) * log(1 - w[ind[i]]) + z[i] * ( log(w[ind[i]]) + logNormal[i] )
Lik[i] <- exp(logLik[i])

# Use the ones trick
p[i] <- Lik[i] / C
ones[i] ~ dbern(p[i])

#y1[i] ~ dpois(x1[ind[i]])

}



# PRIORS
beta0 ~ dnorm(0, 0.0001)
beta1 ~ dnorm(0, 0.0001)
beta2 ~ dnorm(0, 0.0001)
beta3 ~ dnorm(0, 0.0001)
beta4 ~ dnorm(0, 0.0001)
#beta5 ~ dnorm(0, 0.0001)


#gamma0 ~ dnorm(0, 0.0001)
#gamma1 ~ dnorm(0, 0.0001)
#gamma2 ~ dnorm(0, 0.0001)
#gamma3 ~ dnorm(0, 0.0001)
#gamma4 ~ dnorm(0, 0.0001)

#tau2x ~ dgamma(1, 1)
tau2y ~ dgamma(1, 1)

theta0 ~ dnorm(0, 0.0001)
theta1 ~ dnorm(0, 0.0001)
theta2 ~ dnorm(0, 0.0001)

alpha0 ~ dnorm(0, 0.0001)
#alpha1 ~ dnorm(0, 0.0001)
alpha2 ~ dnorm(0, 0.0001)

#shape ~ dgamma(1,1)
shape2 ~ dgamma(1,1)

}"


z1 <- rep(1,nrow(data$y2))
z2 <- rep(1,nrow(data$y2))
z1[data$y2[,1]==0] <- 0
z2[data$y2[,2]==0] <- 0
z <- c(z1,z2)

dattp <- list(y1=c(data$y1[,1],data$y1[,2]),y2=c(data$y2[,1],data$y2[,2]),z=z,age=Za[,2],gender=Za[,3],bmi=Za[,4],smoke=Za[,5],
              N=nrow(data$y2),ones=rep(1,2*nrow(data$y2)),ind=rep(1:nrow(data$y2),2))

mtp = jags.model(textConnection(modeltp), dattp,n.adapt=1000,n.chains=3)
rtp = coda.samples(mtp, c("shape","beta0","beta1","beta2","beta3","beta4","tau2y","gamma0","gamma1","gamma2","gamma3","gamma4","shape2","theta0","theta1"), n.iter=1000)

