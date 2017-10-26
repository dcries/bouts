#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



// Description: random number generator for K dimensional multivariate normal
// n: number of random generations
// mu: mean vector of length K
// sigma: covariance matrix of dim KxK
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}


// Description: random number generator for inverse wishart
// n: number of random generations
// v: degrees of freedom
// S: scale matrix
arma::cube rinvwish(int n, int v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
    arma::mat A(p,p, arma::fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df)); 
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}


const double log2pi = std::log(2.0 * M_PI);

// Description: density of multivariate normal
// x: matrix of values you want density for, each row is one observation
// mean: mean vector for observations
// sigma: covariance matrix for all observations
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// Description: density of multivariate normal
// x: vector of values you want density for
// mean: mean vector for observations
// sigma: covariance matrix for all observations
double dmvnrm_arma(arma::vec x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  int xdim = x.size();
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  

  arma::vec z = rooti * arma::trans( x.t() - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// Description: density of generalized poisson
// x: value to calculate density for
// mu: mu parameter of Gen Poisson
// lambda: lambda parameter of Gen Poisson
// logd: log the density
// [[Rcpp::export]]
double dgenpois(double x, double mu, double lambda, bool logd=true){
  double out;
  out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
  
  if(!logd){
    out=exp(out);
  }
  return out;
}

// Description: density of generalized poisson
// x: value to calculate density for
// mu: vector of mu parameters of Gen Poisson
// lambda: lambda parameter of Gen Poisson
// logd: log the density
arma::vec dgenpois(double x, arma::vec mu, double lambda, bool logd=true){
  int n = mu.size();
  arma::vec out(n);
  out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
  if(!logd){
    out=exp(out);
  }
  return out;
}

// Description:returns vector where element i equals zero when x[i]=0, one otherwise
//x: numeric vector
arma::vec check0(arma::vec x){
  int n = x.size();
  arma::vec out = arma::ones(n);
  for(int i=0;i<n;i++){
    if(x[i]==0){
      out[i] = 0;
    }
  }
  return out; 
}

// Description: subsets x for values of ind=1, x and ind equal length
// x: numeric vector
// ind: vector of ones and zeros of equal length as x
arma::vec subset(arma::vec x, arma::vec ind){
  int n = x.size();
  int m = arma::accu(ind);
  if(x.size() != ind.size()){
    stop("error:check subset function!");
  }
  arma::vec out(m);
  int count = 0;
  for(int i=0;i<n;i++){
    if(ind[i]==1){
      out[count] = x[i];
      count++;
    }
  }
  return out;
}

// Description: subsets x for values of ind=1, x and ind equal length
// x: numeric vector
// ind: vector of ones and zeros of equal length as x
arma::mat subset2(arma::mat x, arma::vec ind){
  int n = x.n_rows;
  int n2 = ind.size();
  int m = arma::accu(ind);
  int k = x.n_cols;
  if(n != n2){
    stop("error:check subset (mat) function!");
  }
  arma::mat out(m,k);
  int count = 0;
  for(int i=0;i<n;i++){
    if(ind[i]==1){
      out.row(count) = x.row(i);
      count++;
    }
  }
  return out;
}

// Description: calculates value of pi (P(Y_1ij > 0|i)) for given x1 and lambda
// x1: vector of mean parameters
// lambda: lambda parameter for Gen Poisson
arma::vec calc_p(arma::vec x1, double lambda){
  int n = x1.size();
  arma::vec out = arma::zeros(n);
  out = 1.0-dgenpois(0,x1,lambda,false);
  return out;
}

// Description: calculates value of pi (P(Y_1ij > 0|i)) for given x1 and lambda
// x1: vector of mean parameters
// lambda: lambda parameter for Gen Poisson
double calc_p(double x1, double lambda){
  double out = 0.0;
  out = 1.0-dgenpois(0,x1,lambda,false);
  return out;
}

// Description:samples beta (regression parameters for average excess MET-mins) from full conditional
// Zbp: model matrix with covariates for those with positive Y_2ij
// mu0b: prior mean for beta
// V0b: prior covariance matrix for beta
// sigma2: variance of log(Y_2ij)
// xp: positive values of Y_2ij
// bp: values of b_2i corresponding to those with positive values of Y_2ij
arma::vec sample_beta(arma::mat Zbp, arma::vec mu0b, arma::mat V0b,
                        double sigma2, arma::vec xp, arma::vec bp){
  int k = Zbp.n_cols;
  arma::vec out(k);
  arma::mat Vb = arma::inv(arma::inv(V0b)+Zbp.t()*Zbp/sigma2);
  arma::vec mub = Vb*(arma::inv(V0b)*mu0b + Zbp.t()*(log(xp)-bp)/sigma2);
  ////std::cout << "Vb=" << sigma2 <<"\n";
  out = (mvrnormArma(1,mub,Vb).row(0).t());
  
  return out;
}


// Description: log likelihood + prior for full conditional of gamma
// cx: prior mean for gamma
// vx: prior covariance matrix for gamma
// lambda: lambda parameter for Gen Poisson
// mu1: mean vector for Y_1
// y1: Y_1 data, number of bouts
// y2: Y_2 data, average excess MET-minutes
// muy2: mean vector for Y_2, in matrix form, same dim as y2
// p: vector of P(Y_1 > 0)
// sigma2y: variance of log(Y_2ij)
// gam: current gamma parameters 
double log_qgamma(arma::vec cx, arma::mat vx, double lambda, arma::vec mu1, 
                  arma::mat y1, arma::mat y2, arma::mat muy2, arma::vec p, double sigma2y, arma::vec gam){
  double out = 0.0;
  int n = y1.n_rows;
  int k = y1.n_cols;
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out += dgenpois(y1(i,j),mu1[i],lambda,true);
      if(y2(i,j)==0){
        out += log(1-p[i]);
      }
      else{
        out += log(p[i]) + R::dlnorm(y2(i,j),muy2(i,j),sqrt(sigma2y),true);
      }
    }
  }
  out += dmvnrm_arma(gam,cx.t(),vx,true);
  return(out);
}

// Description: samples gamma from full conditional via random walk
// cx: prior mean for gamma
// vx: prior covariance matrix for gamma
// lambda: lambda parameter for Gen Poisson
// mu1: mean vector for Y_1
// y1: Y_1 data, number of bouts
// y2: Y_2 data, average excess MET-minutes
// muy2: mean vector for Y_2, in matrix form, same dim as y2
// p: vector of P(Y_1 > 0)
// sigma2y: variance of log(Y_2ij)
// gam: current gamma parameters 
// Za: model matrix 
// Sigma_1: proposal covariance matrix
// currentb: current values of random effects b
// divisor: just in case acceptance rates too high/low
arma::vec sample_gamma(arma::vec cx, arma::mat vx, double lambda, arma::vec mu1, 
                       arma::mat y1, arma::mat y2, arma::mat muy2, arma::mat Za, 
                       arma::vec gamma, arma::vec p, double sigma2y,
                       arma::mat Sigma_1, arma::mat currentb,double divisor=1.0){
  int k = cx.size();
  int n = y1.n_rows;
  double lacceptprob;
  arma::vec out = gamma;
  arma::vec mu1prop(n);
  arma::vec propp(n);
  arma::mat proposalm;
  arma::vec proposal(k);
  proposalm = mvrnormArma(1,gamma,pow(2.4,2)*Sigma_1/(k*divisor));
  
  for(int i=0;i<k;i++){
    proposal[i] = proposalm(0,i);
  }
  
  mu1prop = exp(Za*proposal+currentb.col(0));
  propp = calc_p(mu1prop,lambda);
  
  lacceptprob = log_qgamma(cx,vx,lambda,mu1prop,y1,y2,muy2,propp,sigma2y,proposal)-
    log_qgamma(cx,vx,lambda,mu1,y1,y2,muy2,p,sigma2y,gamma);

  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return(out);
}



// Description: calculates log(mu) for Y_2 
// Za: model matrix for Y_2
// beta: current values of regression coefficients for Y_2
// currentb: current values of b_2, person random effects
arma::vec calc_lmu(arma::mat Za, arma::vec beta, arma::vec currentb){
  int n = Za.n_rows;
  int k = Za.n_cols;
  arma::vec out = arma::zeros(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out[i] += Za(i,j)*beta[j];
    }
    out[i] += currentb[i];
  }
  return out;
}

// Description: samples sigma_y^2 from full conditional
// a0: prior shape parameter for sigma2y
// b0: prior rate for sigma2y
// y1: first column of Y_2, positive values only
// y2: second column of Y_2, positive values only
// Z1: model matrix for above y1, corresponding to positive values only
// Z2: model matrix for above y2, corresponding to positive values only
// beta: current values of regression coefficients beta for Y_2
// bp: current random effects b_2 corresponding to positive values of Y_2
double sample_sigma2(double a0, double b0, arma::vec y1, arma::vec y2, 
                     arma::mat Z1, arma::mat Z2, arma::vec beta, arma::vec bp){
  int n1 = y1.size();
  int n2 = y2.size();
  double out;
  double a = (y1.size()+y2.size())/2.0 + a0;
  double b = b0;
  b += 0.5*arma::accu(pow(log(y1)-Z1*beta-bp.subvec(0,n1-1),2.0));
  b += 0.5*arma::accu(pow(log(y2)-Z2*beta-bp.subvec(n1,n1+n2-1),2.0));
  
  out = 1/R::rgamma(a,1/b);
  return(out);
}  



// Description: log likelihood + prior for lambda
// y1: matrix of observed Y_1
// y2: matrix of observed Y_2
// x1: mean vector for Y_1
// muy2: mean matrix for Y_2
// p: vector of P(Y_1 >0)
// sigma2y: variance of log(Y_2)
// lambda: lambda parameter for Gen Poisson
// a_lambda: Prior for lambda, Beta(a,b)
// b_lambda: Prior for lambda
double log_ql(arma::mat y1, arma::mat y2, arma::vec x1, arma::mat muy2,
              arma::vec p,double sigma2y,double lambda, double a_lambda, double b_lambda){
  double out = 0.0;
  int k = y1.n_cols;
  int n = y1.n_rows;
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out += dgenpois(y1(i,j),x1[i],lambda);
      if(y2(i,j)==0){
        out += log(1-p[i]);
      }
      else{
        out += log(p[i]) + R::dlnorm(y2(i,j),muy2(i,j),sqrt(sigma2y),true);
      }
    }
    
  }
  out += (a_lambda-1)*log(lambda) + (1-b_lambda)*log(1-lambda);
  return out;
}

// Description: sample lambda from full conditional via independence metropolis
// y1: matrix of observed Y_1
// y2: matrix of observed Y_2
// x1: mean vector for Y_1
// muy2: mean matrix for Y_2
// p: vector of P(Y_1 >0)
// sigma2y: variance of log(Y_2)
// lambda: lambda parameter for Gen Poisson
// a_lambda: Prior for lambda, Beta(a,b)
// b_lambda: Prior for lambda
// proplambda1: value for proposal of lambda, Beta(a,b)
// proplambda2: value for proposal of lambda, Beta(a,b)
double sample_lambda(arma::mat y1, arma::mat y2, arma::vec x1, arma::mat muy2,
                     arma::vec p, double sigma2y, double lambda, double a_lambda, double b_lambda,
                     double proplambda1, double proplambda2){
  double out = lambda;
  double proposal = R::rbeta(proplambda1,proplambda2);
  double lacceptprob;
  int n = x1.size();
  arma::vec pprop(n);
  pprop = calc_p(x1,proposal);
  
  lacceptprob = log_ql(y1,y2,x1,muy2,pprop,sigma2y,proposal,a_lambda,b_lambda) - log_ql(y1,y2,x1,muy2,p,sigma2y,lambda,a_lambda,b_lambda) -
    R::dbeta(proposal,proplambda1,proplambda2,true) + R::dbeta(lambda,proplambda1,proplambda2,true);
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return out;
}

// Description: samples Sigma_b from full conditional
// d0: prior for degrees of freedom
// D0: prior scale matrix
// currentb: matrix of current b values, person random effects, first column b_1, second b_2
arma::mat sample_Sigmab(int d0, arma::mat D0, arma::mat currentb){
  int n = currentb.n_rows;
  int df = d0 + n;
  int k = currentb.n_cols;
  arma::mat scale = D0 + currentb.t()*currentb;
  arma::mat out(k,k);
  out = rinvwish(1,df,scale);
  return out;
}

// Description: log likelihood for terms relating to random effects b
// y1: Y_1 values for one individual
// y2: Y_2 values for one individual
// muy: mean parameter for Y_1
// muy2: mean vector parameter for Y_2
// lambdau: lambda parameter for Gen Poiss
// sigma2y: variance of log(Y_2)
// currentb: value of b_1 and b_2 for individual
// Sigmab: covariance matrix for random effects b
double log_qb(arma::rowvec y1, arma::rowvec y2, double muy, arma::vec muy2, double lambdau,
              double sigma2y, arma::vec currentb, arma::mat Sigmab){
  int n = y1.size();
  int kp = currentb.size();
  double ll=0.0;
  double p0;
  arma::rowvec zs = arma::zeros(kp).t();
  
  
  p0 = 1-dgenpois(0,muy,lambdau,false);
  for(int i=0;i<n;i++){
    ll += dgenpois(y1[i],muy,lambdau,true); 
    if(y2[i]==0){
      ll += log(1-p0);
    }
    else{
      ll += log(p0) + R::dlnorm(y2[i],muy2[i],sqrt(sigma2y),true);
    }
  }
  ll += arma::accu(dmvnrm_arma(currentb,zs,Sigmab,true));
  return ll;
}



// samples person random effects b from full conditional via random walk
// y1: Y_1 values for one individual
// y2: Y_2 values for one individual
// muy1: mean vector for Y_1
// muy2: mean matrix parameter for Y_2
// lambdau: lambda parameter for Gen Poiss
// sigma2y: variance of log(Y_2)
// currentb: values of b_1 (first column) and b_2 (second column) for individual
// Sigmab: covariance matrix for random effects b
// gamma: current gamma regression coefficients for Y_1 regression
// betay: current beta regression coefficients for Y_2 regression
// Zg: full model matrix for all individuals
// Zby1: model matrix for those with positive Y_1i1
// Zby2: model matrix for those with positvie Y_1i2
// Sigma_1: proposal covariance matrix 
arma::mat sample_b(arma::mat y1, arma::mat y2,arma::vec muy1,arma::mat muy2,arma::mat currentb, 
                   arma::mat Sigmab, arma::vec gamma, arma::vec betay,
                   arma::mat Zg, arma::mat Zby1, arma::mat Zby2, double lambda, double sigma2y, 
                   arma::cube Sigma_1){

  int n = y1.n_rows;
  int kp = currentb.n_cols;
  int kzg = Zg.n_cols;
  int kzb = Zby1.n_cols;
  double lacceptprob;
  arma::mat out = currentb;
  arma::mat propb2(1,kp);
  arma::vec propb(kp);
  double muy1prop;
  arma::vec muy2prop(2);
  
  for(int i=0;i<n;i++){
    muy1prop = 0.0;
    muy2prop = arma::zeros(2);
    propb2 = mvrnormArma(1,currentb.row(i).t(),pow(2.4,2)*Sigma_1.slice(i)/kp);

    for(int j=0;j<kp;j++){
      propb[j] = propb2(0,j);
    }
    
    for(int j=0;j<kzg;j++){
      muy1prop += Zg(i,j)*gamma[j];
    }
    muy1prop += propb[0];
    
    for(int j=0;j<kzb;j++){
      muy2prop[0] += Zby1(i,j)*betay[j];
      muy2prop[1] += Zby2(i,j)*betay[j];
      
    }
    muy2prop[0] += propb[1];
    muy2prop[1] += propb[1];
    
    lacceptprob = log_qb(y1.row(i),y2.row(i),exp(muy1prop),(muy2prop),lambda,sigma2y,propb,Sigmab) - 
      log_qb(y1.row(i),y2.row(i),muy1[i],muy2.row(i).t(),lambda,sigma2y,currentb.row(i).t(),Sigmab);

        if(log(R::runif(0,1)) < lacceptprob){
      out.row(i) = propb.t();
    }
    
  }
  
  return out;
}

// Description: Run the MCMC for our full model
//data is a list with:
//    y1- a nx2 matrix, n individuals with 2 observations each containing values for number of bouts
//    y2- a nx2 matrix, n individuals with 2 observations each containing values for average excess MET-minutes
//    Za- a nx8 matrix, n individuals with values of covariates
// init is a list with starting values, as well as values for tuning parameters
//    needs: currentb, currentbetay, currentlambda, currentgamma, currentsigma2b, currentSigmab, gammatune, btune, proplambda1,proplambda2
// prior is a list with values for priors
//    needs: m_beta,V_beta,m_gamma,V_gamma,d0,D0,a_sigma2y,b_sigma2y,a_lambda,b_lambda
// nreps is how many MCMC iterations to run
// burn is the number of draws for burn-in
// thin is how many iterations to thin
// [[Rcpp::export]]
List mcmc_2part(     List data, 
                     List init, 
                     List priors, 
                     const int nreps, 
                     const int burn=1000,
                     const int thin=1){
  

  
  //data
  arma::mat Za                = as<arma::mat>(data["Za"]);
  arma::mat y1                = as<arma::mat>(data["y1"]);
  arma::mat y2                = as<arma::mat>(data["y2"]);
  
  int n                       = Za.n_rows;
  int na                      = Za.n_cols;
  arma::vec ybar1             = mean(y1,1);
  arma::vec ybar2             = mean(y2,1);

  
  //starting values
  arma::mat currentb                      = as<arma::mat>(init["currentb"]);
  arma::mat currentSigmab                 = as<arma::mat>(init["currentSigmab"]);
  arma::vec currentbetay                  = as<arma::vec>(init["currentbetay"]);
  arma::vec currentgamma                  = as<arma::vec>(init["currentgamma"]);
  arma::vec currentsigma2b(2);
  currentsigma2b[0]= currentsigma2b[1]    = as<double>(init["currentsigma2b"]);
  double currentsigma2y                   = as<double>(init["currentsigma2y"]);
  double currentlambda                    = as<double>(init["currentlambda"]);

  // proposal values for independence metropolis
  double proplambda1                      = as<double>(init["proplambda1"]);
  double proplambda2                      = as<double>(init["proplambda2"]);

  
  // starting diagonal of covariance matrix for tuning parameters
  arma::vec gammatune                     = as<arma::vec>(init["gammatune"]);
  arma::vec btune                         = as<arma::vec>(init["btune"]);
  
  // storage for current values of mean functions for y1 and y2
  arma::vec currentp(n);
  arma::vec currentlmuy1(n);
  arma::vec currentlmuy2_1(n);
  arma::vec currentlmuy2_2(n);

  
  //priors
  arma::vec m_beta              = as<arma::vec>(priors["m_beta"]);
  arma::vec m_gamma             = as<arma::vec>(priors["m_gamma"]);
  arma::mat V_beta              = as<arma::mat>(priors["V_beta"]);
  arma::mat V_gamma             = as<arma::mat>(priors["V_gamma"]);
  arma::mat D0                  = as<arma::mat>(priors["D0"]);
  double a_sigma2y              = as<double>(priors["a_sigma2y"]);
  double b_sigma2y              = as<double>(priors["b_sigma2y"]);
  double a_lambda               = as<double>(priors["a_lambda"]);
  double b_lambda               = as<double>(priors["b_lambda"]);
  double d0                     = as<double>(priors["d0"]);

  // number of random effects per person
  int nbre = currentb.n_cols;
  // number of MCMC iterations to store
  int keep = (nreps-burn)/thin;
  
  //storage
  arma::mat betay(keep,m_beta.size());
  arma::mat gamma(keep,na);
  arma::vec sigma2y(keep);
  arma::vec lambda(keep);
  arma::mat b1(keep,n);
  arma::mat b2(keep,n);
  arma::mat sigma2b(keep,nbre);
  arma::vec corrb(keep);

  // holds values of burnin to adapt covariance matrix of random walk
  arma::mat gammaburn(burn,na);
  arma::mat b1burn(burn,n);
  arma::mat b2burn(burn,n);
  
  
  arma::mat Za1positive;
  arma::mat Za2positive;
  arma::mat Zapositive;
  
  arma::mat gamma_var = arma::zeros(na,na);
  gamma_var.diag() = gammatune;
  
  arma::cube b_var = arma::zeros(nbre,nbre,n);
  for(int i=0;i<n;i++){
    b_var.slice(i).diag() = btune;
  }

  arma::vec indy1 = check0(y2.col(0));
  arma::vec indy2 = check0(y2.col(1));
  arma::vec indybar = check0(ybar2);
  
  arma::vec y2sub1 = subset(y2.col(0),indy1); // values of y2 s.t. first obs >0
  arma::vec y2sub2 = subset(y2.col(1),indy2); // values of y2 s.t. second obs >0
  arma::vec y2sub = arma::join_cols(y2sub1,y2sub2);

  Za1positive = subset2(Za,indy1); //rows of Za s.t. first observation of y1>0
  Za2positive = subset2(Za,indy2);//rows of Za s.t. second observation of y1>0
  Zapositive = arma::join_cols(Za1positive,Za2positive);
  
  arma::mat b;

  // vector of b s.t. y1 is >0
  arma::vec b_positivey1 = arma::join_cols(subset(currentb.col(1),indy1),subset(currentb.col(1),indy2));
  currentlmuy2_1 = calc_lmu(Za,currentbetay,currentb.col(1));
  currentlmuy2_2 = calc_lmu(Za,currentbetay,currentb.col(1));
  
  // sets index for storage
  int ind = 0;
  
  for(int i=0;i<nreps;i++){

    currentbetay = sample_beta(Zapositive,m_beta,V_beta,currentsigma2y,y2sub,b_positivey1);
    
    currentlmuy2_1 = calc_lmu(Za,currentbetay,currentb.col(1));
    currentlmuy2_2 = calc_lmu(Za,currentbetay,currentb.col(1));
    
    currentgamma = sample_gamma(m_gamma,V_gamma,currentlambda,exp(currentlmuy1),
                                y1,y2,arma::join_rows(currentlmuy2_1,currentlmuy2_2),
                                Za,currentgamma,currentp,currentsigma2y,gamma_var,currentb,4.0);

    currentlmuy1 = calc_lmu(Za,currentgamma,currentb.col(0));
    
    currentlambda = sample_lambda(y1,y2,exp(currentlmuy1),arma::join_rows(currentlmuy2_1,currentlmuy2_2),
                                  currentp,currentsigma2y,currentlambda,a_lambda,b_lambda,proplambda1,proplambda2);
    
    currentp = calc_p(exp(currentlmuy1),currentlambda);

    currentSigmab = sample_Sigmab(d0,D0,currentb);

    
    if((i>99)&&(i<burn)&&(i%20==0)){
      gamma_var = cov(gammaburn.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(gamma_var(j,j)==0){
          std::cout << "gamma proposal covariance matrix has variance 0\n";
        }
      }
    }
    
    currentb = sample_b(y1,y2,exp(currentlmuy1),arma::join_rows(currentlmuy2_1,currentlmuy2_2),
                        currentb,currentSigmab,
                        currentgamma,currentbetay,Za,Za,Za,currentlambda,
                        currentsigma2y,b_var);
    
    currentlmuy1 = calc_lmu(Za,currentgamma,currentb.col(0));
    
    currentp = calc_p(exp(currentlmuy1),currentlambda);
    
    currentlmuy2_1 = calc_lmu(Za,currentbetay,currentb.col(1));
    
    currentlmuy2_2 = calc_lmu(Za,currentbetay,currentb.col(1));
    
    if((i>99)&&(i<burn)&&(i%20==0)){
      for(int j=0;j<n;j++){
        b = arma::join_rows(b1burn.col(j),b2burn.col(j));
        b_var.slice(j) = cov(b.rows(0,i-1));
      }
    }
    
    
    b_positivey1 = arma::join_cols(subset(currentb.col(1),indy1),subset(currentb.col(1),indy2));

    currentsigma2y = sample_sigma2(a_sigma2y,b_sigma2y,y2sub1,y2sub2,Za1positive,Za2positive,currentbetay,b_positivey1);
    
    // store values during burnin so we can update proposal covariance matricies
    if(i < burn){
      gammaburn.row(i) = currentgamma.t();    
      b1burn.row(i) = currentb.col(0).t();
      b2burn.row(i) = currentb.col(1).t();
    }
    // store data after burnin and only on appropriate iterations as determined by thinning
    if((i >= burn) && (i%thin==0) ){
      betay.row(ind) = currentbetay.t();
      gamma.row(ind) = currentgamma.t();    
      sigma2y[ind] = currentsigma2y;
      lambda[ind]      = currentlambda;
      b1.row(ind) = currentb.col(0).t();
      b2.row(ind) = currentb.col(1).t();
      sigma2b.row(ind) = (currentSigmab.diag()).t();
      corrb[ind] = currentSigmab(0,1)/(sqrt(currentSigmab(0,0)*currentSigmab(1,1)));
      ind++;
    }
    
 
    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  } 
  
  return List::create(
    Named("betay") = betay,
    Named("gamma") = gamma,
    Named("sigma2y") = sigma2y,
    Named("lambda")    = lambda,
    Named("b1") = b1,
    Named("b2") = b2,
    Named("sigma2b") = sigma2b,
    Named("corrb") = corrb);
} 




