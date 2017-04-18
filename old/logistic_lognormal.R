library(rjags)
library(rstan)

pams <- read.csv("C:\\Users\\dcries\\github\\pams\\FinalPAMSDataSetNew.csv")
pams <- pams[-which(is.na(pams$ModPAR)|is.na(pams$Age)|is.na(pams$BMI)|is.na(pams$Smoker)),]

ind <- rep(0,nrow(pams))
ind[pams$VigPAR > 0] <- 1
pm1 <- glm(ind~Age+Weight+Gender+Smoker+BMI,family="binomial",data=(pams))


n <- 1000
beta1 <- c(0.1,-1.2,4)
beta2 <- c(.6,-0.8,2)
sigmae <- .10
x1 <- matrix(c(runif(n,0,5),runif(n,0,5),rbinom(n,1,0.8)),ncol=3,byrow=FALSE)
x2 <- matrix(c(runif(n,0,5),runif(n,0,5),rbinom(n,1,0.8)),ncol=3,byrow=FALSE)

p <- 1/(1+exp(-x1%*%beta1))
plot(p)
logy <- rnorm(n,x2%*%beta2,sigmae)
y <- exp(logy)
y[x2[,3]==0] <- 0
plot(y)

z <- rep(1,length(y))
z[y==0] <- 0

m1 <- glm(z~x1+0,family="binomial")


models <- '
data {
   int<lower=0> N ;
   int p;
   vector<lower=0>[N] y ; // data
   //vector[N] z; //0 or 1
   matrix[N,p] x;
   vector[p] m;     //prior mean vector for betas
   matrix[p,p] S;   //prior covariance for betas
   real nu;         //prior for degrees of freedom? matrix for random effects
   matrix[2,2] lambda; //prior for cov matrix for random effects
   vector[2] zero; //c(0,0) for prior mean of random effects
}
parameters {
  real<lower=0> sigmae ; // shape
  vector[p] beta1;
  vector[p] beta2;
  vector[2] b;
  cov_matrix[2] Sb;

}
transformed parameters{
  vector<lower=0, upper=1>[N] pi ; // prob.
  vector<lower=0>[N] mu ; // prob.
  vector[N] temp1;
  vector[N] temp2;
  real<lower=0> sb1;
  real<lower=0> sb2;
  real<lower=-1,upper=1> corb;

  sb1 <- pow(Sb[1,1],0.5);
  sb2 <- pow(Sb[2,2],0.5);
  corb <- Sb[1,2]/(sb1*sb2);

  temp1 <- x*beta1+b[1];
  temp2 <- x*beta2+b[2];
  for(i in 1:N){
    //pi[i] <- 1/(1+pow(e(),-temp1[i]));
    pi[i] <- inv_logit(temp1[i]);
    //mu[i] <- pow(e(),temp2[i]);
    mu[i] <- exp(temp2[i]);
  }

}
model {
//  for (i in 1:N) {
//    if ( y[i] == 0 )
//      increment_log_prob(bernoulli_log(0,pi[i]));
//    else
//      increment_log_prob(bernoulli_log(1,pi[i]) + lognormal_log(y[i],mu[i],sigmae));
//  }
    for(n in 1:N){
        increment_log_prob( if_else(y[n] == 0, log1m(pi[n]), log(pi[n]) + lognormal_log(y[n], mu[n], sigmae)) );
    }

  b ~ multi_normal(zero,Sb);//random effect
  Sb ~ inv_wishart(nu,lambda);
  sigmae ~ cauchy(0,1);
  beta1 ~ multi_normal(m,S);
  beta2 ~ multi_normal(m,S);

}'

datalist <- list(y=y,
            x=x1,
            N=length(y),
            p=ncol(x1),
            m=rep(0,ncol(x1)),
            S=10*diag(ncol(x1)),
            nu=3,
            lambda=diag(2),
            zero=rep(0,2))

x1 <- with(pams,model.matrix(~Gender))
x2 <- scale(cbind(pams$Age,pams$BMI),center = T,scale = T)
x <- cbind(x1,x2)
datalist <- list(y=pams$ModPAR,
            x=x,
                  N=nrow(pams),
                  p=ncol(x),
                  m=rep(0,ncol(x)),
                  S=10*diag(ncol(x)),
                  nu=3,
                  lambda=diag(2),
                  zero=rep(0,2))

model2 = stan_model(model_code=models)
results2 = sampling(model2,data=datalist,pars=c('beta1','beta2','sigmae','sb1','sb2','corb'),chains=1)#,init=list(list(beta1=c(0,0,0),beta2=c(0,0,0),sigmae=1)))

