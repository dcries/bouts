#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "rtn1.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}


const double log2pi = std::log(2.0 * M_PI);

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
    //std::cout << x.row(i) << "\n" << mean;
    
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}
// // [[Rcpp::export]]
double dmvnrm_arma(arma::vec x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //int n = x.n_rows;
  
  int xdim = x.size();
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  //for (int i=0; i < n; i++) {
  //std::cout << x.row(i) << "\n" << mean;
  
  arma::vec z = rooti * arma::trans( x.t() - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;
  //}
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

double dgenpois(double x, double mu, double lambda){
  double out;
  out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
  return out;
}

//does currentx > 0
// [[Rcpp::export]]
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

// // [[Rcpp::export]]
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

// [[Rcpp::export]]
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

//// [[Rcpp::export]]
arma::vec sample_u(arma::mat Za, arma::vec x, arma::vec alpha){
  int n = Za.n_rows;
  int k = Za.n_cols;
  arma::vec out(n);
  //double check;
  double mu;
  for(int i=0;i<n;i++){
    mu = 0.0;
    for(int j=0;j<k;j++){
      mu += Za(i,j)*alpha[j];
    }
    //check = R::rnorm(mu+b(i,0),1); 
    //std::cout << i << "\n";
    if(x[i]==0){
      //out[i] = -std::abs(check);
      out[i] = -rtn1(mu,1,0,999999);
    }
    else{
      //out[i] = std::abs(check);
      out[i] = rtn1(mu,1,0,999999);
    }
  }
  return out;
}

//// [[Rcpp::export]]
arma::vec sample_alpha(arma::mat Za,  arma::vec u, 
                       arma::vec mu0a, arma::mat V0a){
  int k = Za.n_cols;
  arma::vec out(k);
  arma::mat Va = arma::inv(arma::inv(V0a)+Za.t()*Za);
  arma::vec mua = Va*(arma::inv(V0a)*mu0a + Za.t()*(u));
  
  out = (mvrnormArma(1,mua,Va).row(0).t());
  
  
  return out;
}

//// [[Rcpp::export]]
arma::vec calc_p(arma::mat Za, arma::vec alpha){
  int n = Za.n_rows;
  int k = Za.n_cols;
  arma::vec out(n);
  double x;
  for(int i=0;i<n;i++){
    x=0.0;
    for(int j=0;j<k;j++){
      x += Za(i,j)*alpha[j];
    }
    out[i] = R::pnorm5(x,0,1,1,0);
  }
  return out;
}

//// [[Rcpp::export]]
arma::vec sample_beta_2(arma::mat Zbp, arma::vec mu0b, arma::mat V0b,
                      double sigma2, arma::vec xp){
  int k = Zbp.n_cols;
  arma::vec out(k);
  arma::mat Vb = arma::inv(arma::inv(V0b)+Zbp.t()*Zbp/sigma2);
  arma::vec mub = Vb*(arma::inv(V0b)*mu0b + Zbp.t()*(log(xp))/sigma2);
  //std::cout << "Vb=" << sigma2 <<"\n";
  out = (mvrnormArma(1,mub,Vb).row(0).t());
  
  return out;
}


// [[Rcpp::export]]
double log_qgamma(arma::vec cx, arma::mat vx, double eta, arma::vec mu1, arma::vec x1, arma::vec gam){
  double out = 0.0;
  int n = x1.size();
  for(int i=0;i<n;i++){
    out += R::dgamma(x1[i],eta,1/(eta/mu1[i]),true);
  }
  out += dmvnrm_arma(gam,cx.t(),vx,true);
  return(out);
}

// [[Rcpp::export]]
arma::vec sample_gamma(arma::vec cx, arma::mat vx, double eta, arma::vec mu1, 
                       arma::vec x1, arma::mat Za, arma::vec gamma, 
                       arma::mat Sigma_1, double divisor=1.0){
  int k = cx.size();
  int n = x1.size();
  double lacceptprob;
  arma::vec out = gamma;
  arma::vec mu1prop(n);
  arma::mat proposalm;
  arma::vec proposal(k);
  // std::cout << "before mvrnorm \n";
  // std::cout << "gamma= " << gamma << "\n";
  // std::cout << "Sigma= " << Sigma_1 << "\n";
  proposalm = mvrnormArma(1,gamma,pow(2.4,2)*Sigma_1/(k*divisor));
  
  //std::cout << "after mvrnorm\n";
  for(int i=0;i<k;i++){
    proposal[i] = proposalm(0,i);
  }

  mu1prop = exp(Za*proposal);

  lacceptprob = log_qgamma(cx,vx,eta,mu1prop,x1,proposal)-log_qgamma(cx,vx,eta,mu1,x1,gamma);
  //std::cout << "after dvrnorm\n";
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return(out);
}


//// [[Rcpp::export]]
arma::vec calc_lmu(arma::mat Zb, arma::vec beta){
  int n = Zb.n_rows;
  int k = Zb.n_cols;
  arma::vec out = arma::zeros(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out[i] += Zb(i,j)*beta[j];
    }
  }
  return out;
}

//// [[Rcpp::export]]
double sample_sigma2(double a0, double b0, arma::vec xp, arma::mat Zbp,
                     arma::vec beta){
  double out;
  double a = xp.size()/2.0 + a0;
  double b = b0 + 0.5*arma::accu(pow(log(xp)-Zbp*beta,2.0));
  out = 1/R::rgamma(a,1/b);
  return(out);
}  



double sample_sigma2(double a0, double b0, arma::vec y1, arma::vec y2, 
                     arma::mat Z1, arma::mat Z2, arma::vec beta){

  double out;
  double a = (y1.size()+y2.size())/2.0 + a0;
  double b = b0;
  b += 0.5*arma::accu(pow(log(y1)-Z1*beta,2.0));
  b += 0.5*arma::accu(pow(log(y2)-Z2*beta,2.0));
                                              
  out = 1/R::rgamma(a,1/b);
  return(out);
}  

double log_ex(double ae, double be, double eta, arma::vec x1, arma::vec mu1){
  double out = (ae-1)*log(eta) - eta*be;
  int n = x1.size();
  for(int i=0;i<n;i++){
    out += R::dgamma(x1[i],eta,1/(eta/mu1[i]),true);
  }
  return(out);
}

// [[Rcpp::export]]
double sample_eta(double ae, double be, double currenteta, arma::vec x1, arma::vec mu1, double propa, double propb){
  double proposal;
  double lacceptprob;
  double out = currenteta;

  proposal = R::rgamma(propa,1/propb);
  lacceptprob = log_ex(ae,be,proposal,x1,mu1) - log_ex(ae,be,currenteta,x1,mu1) -
    R::dgamma(proposal,propa,1/propb,true) + R::dgamma(currenteta,propa,1/propb,true);
  
  if(log(R::runif(0,1)) < lacceptprob){
    out = proposal;
  }
  return(out);
}


// [[Rcpp::export]]
double log_qx1(double x1, double x2, arma::rowvec y1, arma::rowvec y2, double p, double muy2,
               double mux1, double mux2,  double sigma2x, double sigma2y, double eta, double lambda){
  double ll1;
  double ll2;
  double ll3=0.0;
  double ll4=0.0;
  
  int k = y2.size();
  
  //ll1 = R::dlnorm(x2,mux2,sqrt(sigma2x),true);
  ll1 = R::dgamma(x2,sigma2x,1.0/(sigma2x/mux2),true);
  ll2 = R::dgamma(x1,eta,1.0/(eta/mux1),true);
  
  for(int i=0;i<k;i++){
    //ll4 += R::dpois(y1[i],x1,true);
    ll4 += dgenpois(y1[i],x1,lambda);
    if(y2[i]==0){
      ll3 += log(1-p);
    }
    else{
      ll3 += log(p) + R::dlnorm(y2[i],muy2,sqrt(sigma2y),true);
      //ll3 += log(p) + R::dgamma(y2[i],sigma2y,1.0/(sigma2y/muy2),true);
    }
  }
  
  return ll1+ll2+ll3+ll4;
}
// [[Rcpp::export]]
double log_qx2(double x, arma::rowvec y, double p, double muy2,
               double mux2,  double sigma2x, double sigma2y){
  
  if(x <= 0){
    return -1*arma::datum::inf;
  }
  else{
    double ll1;
    double ll3=0.0;
    
    int k = y.size();
    
    //ll1 = R::dlnorm(x,mux2,sqrt(sigma2x),true);
    ll1 = R::dgamma(x,sigma2x,1.0/(sigma2x/mux2),true);
    
    for(int i=0;i<k;i++){
      if(y[i]==0){
        ll3 += log(1-p);
      }
      else{
        ll3 += log(p) + R::dlnorm(y[i],muy2,sqrt(sigma2y),true);
        //ll3 += log(p) + R::dgamma(y[i],sigma2y,1.0/(sigma2y/muy2),true);
      }
    }
    
    return ll1+ll3;
  }
}


// [[Rcpp::export]]
arma::vec sample_x2(arma::vec x1, arma::vec x2, arma::mat y2, arma::vec p, arma::vec mux2, 
                    arma::vec muy2, arma::vec betay, arma::vec alpha, double sigma2x, double sigma2y,
                    double propx2, arma::vec vx2, arma::mat Z){
  int n = x2.size();
  int k = alpha.size();
  double muyprop;
  double pprop;
  double pprop2;
  arma::vec out = x2; 
  
  
  double propx;
  //double propx3;
  
  double lacceptprob;
  
  for(int i=0;i<n;i++){
    //propx = R::rexp(propx2);
    propx = R::rnorm(x2[i],2.4*sqrt(vx2[i]));
    
    // if(propx < 0.0){
    //   propx3 = abs(propx);
    // }
    // else{
    //   propx3 = propx;
    // }

    muyprop = betay[0] + betay[1]*x1[i] + betay[2]*log(propx);
    //muyprop = betay[0]  + betay[1]*propx;
    
    //pprop = R::pnorm5(alpha[0] + alpha[1]*x1[i] + alpha[2]*propx,0,1,1,0);
    pprop2 = 0.0;
    for(int j=0;j<k-1;j++){
      pprop2 += Z(i,j)*alpha[j];
    }
    pprop2 += propx*alpha[k-1];
    pprop = R::pnorm5(pprop2,0,1,1,0);
    
    //pprop = R::pnorm5(alpha[0] + alpha[1]*propx,0,1,1,0);
    
    // if(i==593){
    //   std::cout << "propx = " << propx << " pprop = " << pprop << " muyprop = " <<
    //     muyprop << "x2 = " << x2[i] << " p = " << p[i] << "muy2 = " << muy2[i] << "\n";
    // }
    
    lacceptprob = log_qx2(propx,y2.row(i),pprop,muyprop,mux2[i],sigma2x,sigma2y) -
     log_qx2(x2[i],y2.row(i),p[i],muy2[i],mux2[i],sigma2x,sigma2y);
    
    // lacceptprob = log_qx2(propx,y2.row(i),pprop,muyprop,mux2[i],sigma2x,sigma2y) -
    //   log_qx2(x2[i],y2.row(i),pprop,muyprop,mux2[i],sigma2x,sigma2y);
    
    
    if(log(R::runif(0,1))<lacceptprob){
      out[i]=propx;
    }
  }
  return out;
}

double log_qa(arma::mat y2, arma::vec x2, arma::mat Z,
              arma::vec alpha, double delta,
              arma::vec mu0a, arma::mat V0a){
  double ll=0.0;
  int n = y2.n_rows;
  int p = alpha.size();
  int k = y2.n_cols;
  
  arma::vec est = arma::zeros(n);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<p;j++){
      est[i] += Z(i,j)*alpha[j];
    }
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(y2(i,j)==0){
        ll += log(1-R::pnorm5(est[i],0,1,1,0));
      }
      else{
        ll += R::pnorm5(est[i],0,1,1,1) + R::dgamma(y2(i,j),delta,1.0/(delta/x2[i]),1);
      }
    }
  }
  ll += dmvnrm_arma(alpha,mu0a.t(),V0a,true);
  
  return ll;
}

arma::vec sample_alpha2(arma::mat y2, arma::vec x2, arma::mat Z,
                        arma::vec alpha, double delta, arma::mat sigma_alpha,
                        arma::vec mu0a, arma::mat V0a){
  
  int p = alpha.size();
  double lacceptprob;
  arma::vec out = alpha;
  arma::vec proposal(p);
  proposal = (mvrnormArma(1,alpha,pow(2.4,2)*sigma_alpha/p).row(0)).t();
  
  lacceptprob = log_qa(y2,x2,Z,proposal,delta,mu0a,V0a) -
    log_qa(y2,x2,Z,alpha,delta,mu0a,V0a);
  
  if(log(R::runif(0,1)) < lacceptprob){
    out = proposal;
  }
  
  return out;
}

double log_ed(double ad, double bd, double delta, arma::mat y2, arma::vec x2){
  double out = (ad-1)*log(delta) - delta*bd;
  int n = x2.size();
  int k = y2.n_cols;
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(y2(i,j)>0){
        out += R::dgamma(y2(i,j),delta,1.0/(delta/x2[i]),true);
      }
    }
  }
  return(out);
}

double sample_delta(double ad, double bd, double currentdelta, arma::mat y2, arma::vec x2, double propa, double propb){
  double proposal;
  double lacceptprob;
  double out = currentdelta;
  
  proposal = R::rgamma(propa,1.0/propb);
  lacceptprob = log_ed(ad,bd,proposal,y2,x2) - log_ed(ad,bd,currentdelta,y2,x2) -
    R::dgamma(proposal,propa,1.0/propb,true) + R::dgamma(currentdelta,propa,1.0/propb,true);
  
  if(log(R::runif(0,1)) < lacceptprob){
    out = proposal;
  }
  return(out);
}


//// [[Rcpp::export]]
arma::vec sample_x1(arma::vec x1, arma::vec x2, arma::mat y1, arma::mat y2, 
                    arma::mat z, arma::vec p, arma::vec mux1, arma::vec mux2, 
                    arma::vec muy, arma::vec betay, arma::vec betax, 
                    arma::vec alpha, double sigma2x, double sigma2y,double eta,
                    arma::vec x1propa, arma::vec x1propb, double lambda, arma::mat Zalpha){
  int n = x2.size();
  int k = z.n_cols;
  int kp = alpha.size();
  double muyprop;
  double pprop;
  double pprop2;
  double muxprop;
  double zb;
  arma::vec out = x1; 
  

  double propx;
  double lacceptprob;
  
  for(int i=0;i<n;i++){
    //propx = R::rexp(x1tune[i]);
    propx = R::rgamma(x1propa[i],1/x1propb[i]);
    
    // muyprop = 0.0;
    // for(int j=0;j<k;j++){
    //   muyprop = 
    // }
    zb = 0.0;
    for(int j=0;j<k-1;j++){
      zb += z(i,j)*betax[j];
    }
    

     muyprop = betay[0] + betay[1]*propx + betay[2]*log(x2[i]);
     //pprop = R::pnorm5(alpha[0] + alpha[1]*propx + alpha[2]*x2[i],0,1,1,0);
     
     pprop2 = 0.0;
     for(int j=0;j<kp-2;j++){
       pprop2 += Zalpha(i,j)*alpha[j];
     }
     pprop2 += propx*alpha[kp-2]+Zalpha(i,kp-1)*alpha[kp-1];
     pprop = R::pnorm5(pprop2,0,1,1,0);
     
    //muyprop = betay[0]  + betay[1]*x2[i];
    //pprop = R::pnorm5(alpha[0]  + alpha[1]*x2[i],0,1,1,0);
    
    muxprop = exp(betax[k-1]*propx + zb);

    lacceptprob = log_qx1(propx,x2[i],y1.row(i),y2.row(i),pprop,muyprop,mux1[i],muxprop,sigma2x,sigma2y,eta,lambda) - 
      log_qx1(x1[i],x2[i],y1.row(i),y2.row(i),p[i],muy[i],mux1[i],mux2[i],sigma2x,sigma2y,eta,lambda) -
      //R::dexp(propx,x1tune[i],true) + R::dexp(x1[i],x1tune[i],true);
      R::dgamma(propx,x1propa[i],1/x1propb[i],true) + R::dgamma(x1[i],x1propa[i],1/x1propb[i],true);
    
    if(log(R::runif(0,1))<lacceptprob){
      out[i]=propx;
    }
  }
  return out;
}

double log_ql(arma::mat y1, arma::vec x1, double lambda, double a0l, double b0l){
  double out = 0.0;
  int k = y1.n_cols;
  int n = y1.n_rows;
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out += dgenpois(y1(i,j),x1[i],lambda);
    }
  }
  out += (a0l-1)*log(lambda) + (1-b0l)*log(1-lambda);
  return out;
}

double sample_lambda(arma::mat y1, arma::vec x1, double lambda, double a0l, double b0l,
                     double propl1, double propl2){
  double out = lambda;
  double proposal = R::rbeta(propl1,propl2);
  double lacceptprob;
  
  lacceptprob = log_ql(y1,x1,proposal,a0l,b0l) - log_ql(y1,x1,lambda,a0l,b0l) -
    R::dbeta(proposal,propl1,propl2,true) + R::dbeta(lambda,propl1,propl2,true);
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return out;
}

// [[Rcpp::export]]
List mcmc_2part_1(List data, 
                List init, 
                List priors, 
                const int nreps, 
                const int burn=1000){
  
  //data
  arma::mat Za              = as<arma::mat>(data["Za"]);
  arma::mat Zb              = as<arma::mat>(data["Zb"]);
  arma::mat y1               = as<arma::mat>(data["y1"]);
  arma::mat y2               = as<arma::mat>(data["y2"]);
  
  //..std::cout << "1\n";
  
  int n                     = Za.n_rows;
  //int nr                    = y1.n_cols;
  int na                    = Za.n_cols;
  int nb                    = Zb.n_cols;
  arma::vec ybar1            = mean(y1,1);
  arma::vec ybar2            = mean(y2,1);
  arma::vec intercept       = arma::ones(n);
  
  //std::cout << "2\n";
  
  
  //starting values
  arma::vec currentbetay     = as<arma::vec>(init["currentbetay"]);
  arma::vec currentbetax     = as<arma::vec>(init["currentbetax"]);
  arma::vec currentalpha    = as<arma::vec>(init["currentalpha"]);
  arma::vec currentgamma    = as<arma::vec>(init["currentgamma"]);
  double currentsigma2x     = as<double>(init["currentsigma2x"]);
  double currentsigma2y     = as<double>(init["currentsigma2y"]);
  double currenteta     = as<double>(init["currenteta"]);
  double currenttheta     = as<double>(init["currenttheta"]);
  double currentlambda     = as<double>(init["currentlambda"]);
  double currentdelta     = as<double>(init["currentdelta"]);
  
  
  double propa     = as<double>(init["propa"]);
  double propb     = as<double>(init["propb"]);
  double propax2     = as<double>(init["propax2"]);
  double propbx2     = as<double>(init["propbx2"]);
  double propx2     = as<double>(init["propx2"]);
  double propl1     = as<double>(init["propl1"]);
  double propl2     = as<double>(init["propl2"]);
  double propd1     = as<double>(init["propd1"]);
  double propd2     = as<double>(init["propd2"]);
  
  arma::vec currentx1        = as<arma::vec>(init["currentx1"]);
  arma::vec currentx2        = as<arma::vec>(init["currentx2"]);
  arma::vec currentu        = arma::zeros(n);
  
  arma::vec gammatune            = as<arma::vec>(init["gammatune"]);
  arma::vec betaxtune            = as<arma::vec>(init["betaxtune"]);
  arma::vec alphatune            = as<arma::vec>(init["alphatune"]);
  
  arma::vec vx2             = as<arma::vec>(init["vx2"]);
  //arma::vec x1tune          = as<arma::vec>(init["x1tune"]);
  arma::vec x1propa          = as<arma::vec>(init["x1propa"]);
  arma::vec x1propb          = as<arma::vec>(init["x1propb"]);
  
  arma::vec currentp(n);
  arma::vec currentlmuy(n);
  arma::vec currentlmux1(n);
  arma::vec currentlmux2(n);
  
  //std::cout << "3\n";
  
  
  //priors
  arma::vec mu0y2            = as<arma::vec>(priors["mu0y2"]);
  arma::vec mu0x1            = as<arma::vec>(priors["mu0x1"]);
  arma::vec mu0x2            = as<arma::vec>(priors["mu0x2"]);
  arma::vec mu0a            = as<arma::vec>(priors["mu0a"]);
  arma::mat V0y2             = as<arma::mat>(priors["V0y2"]);
  arma::mat V0x1             = as<arma::mat>(priors["V0x1"]);
  arma::mat V0x2             = as<arma::mat>(priors["V0x2"]);
  arma::mat V0a             = as<arma::mat>(priors["V0a"]);
  double a0eta                 = as<double>(priors["a0eta"]);
  double b0eta                 = as<double>(priors["b0eta"]);
  double a0theta                 = as<double>(priors["a0theta"]);
  double b0theta                 = as<double>(priors["b0theta"]);
  double a0x                = as<double>(priors["a0x"]);
  double b0x                = as<double>(priors["b0x"]);
  double a0y                = as<double>(priors["a0y"]);
  double b0y                = as<double>(priors["b0y"]);
  double a0l                = as<double>(priors["a0l"]);
  double b0l                = as<double>(priors["b0l"]);
  double a0delta                = as<double>(priors["a0delta"]);
  double b0delta                = as<double>(priors["b0delta"]);
  //std::cout << "4\n";
  
  //storage
  arma::mat betay(nreps,mu0y2.size());
  arma::mat betax(nreps,nb+1);
  arma::mat gamma(nreps,nb);
  arma::mat alpha(nreps,mu0a.size());
  arma::vec sigma2x(nreps);
  arma::vec sigma2y(nreps);
  arma::vec eta(nreps);
  arma::vec theta(nreps);
  arma::vec delta(nreps);
  arma::vec lambda(nreps);
  arma::mat latentx1(nreps,n);
  arma::mat latentx2(nreps,n);
  arma::mat muy(nreps,n);
  arma::mat mux1(nreps,n);
  arma::mat mux2(nreps,n);
  arma::mat pg0(nreps,n);
  
  //std::cout << "5\n";
  
  arma::mat x1x2;
  arma::mat x1x2y;
  arma::mat x1x2p1;
  arma::mat x1x2p2;
  arma::mat x1x2p;

  arma::mat gamma_var = arma::zeros(na,na);
  gamma_var.diag() = gammatune;
  
  arma::mat betax_var = arma::zeros(nb+1,nb+1);
  betax_var.diag() = betaxtune;
  
  arma::mat alpha_var = arma::zeros(mu0a.size(),mu0a.size());
  alpha_var.diag() = alphatune;
  
  //std::cout << "5\n";

  arma::vec indy1 = check0(y2.col(0));
  arma::vec indy2 = check0(y2.col(1));
  arma::vec indybar = check0(ybar2);
  
  arma::vec y2sub1 = subset(y2.col(0),indy1);
  arma::vec y2sub2 = subset(y2.col(1),indy2);
  arma::vec y2sub = arma::join_cols(y2sub1,y2sub2);

  //std::cout << "6\n";

  
  
  arma::mat Zbx;
  //std::cout << "2\n";
  //std::cout << "7\n";
  Zbx = arma::join_rows(Zb,currentx1);
  x1x2 = arma::join_rows(intercept,arma::join_rows(currentx1,currentx2));
  
  x1x2y = arma::join_rows(intercept,arma::join_rows(currentx1,log(currentx2)));//x1x2;//arma::join_rows(intercept,currentx2);
  x1x2p1 = subset2(x1x2y,indy1);
  x1x2p2 = subset2(x1x2y,indy2);
  x1x2p = arma::join_cols(x1x2p1,x1x2p2);
  
  
  
  for(int i=0;i<nreps;i++){
    //std::cout << "4\n";
    
    currentu = sample_u(arma::join_cols(x1x2,x1x2), arma::join_cols(y2.col(0),y2.col(1)), currentalpha);
    //std::cout << "5\n";
    
    currentalpha = sample_alpha(arma::join_cols(x1x2,x1x2),currentu,mu0a,V0a);
    //currentalpha = sample_alpha2(y2,currentx2,x1x2,currentalpha,currentdelta,alpha_var,mu0a,V0a);
    currentp = calc_p(x1x2,currentalpha);
    //std::cout << "6\n";
    
    currenteta = sample_eta(a0eta,b0eta,currenteta,currentx1,exp(currentlmux1),propa,propb);
    //std::cout << "6a\n";
    currenttheta = sample_eta(a0theta,b0theta,currenttheta,currentx2,exp(currentlmux2),propax2,propbx2);
    //currentdelta = sample_delta(a0delta,b0delta,currentdelta,y2,currentx2,propd1,propd2);
    
    //////currentbetay = sample_beta_2(x1x2p,mu0y2,V0y2,currentsigma2y,y2sub);
    currentlmuy = calc_lmu(x1x2y,currentbetay);
    
    //std::cout << "6b\n";
    
    //currentbetax = sample_beta_2(Zbx,mu0x2,V0x2,currentsigma2x,currentx2);
    currentbetax = sample_gamma(mu0x2,V0x2,currenttheta,exp(currentlmux2),currentx2,Zbx,currentbetax,betax_var,5.0);
    
    currentlmux2 = calc_lmu(Zbx,currentbetax);
    
    //std::cout << "6c\n";
    
    currentgamma = sample_gamma(mu0x1,V0x1,currenteta,exp(currentlmux1),currentx1,Za,currentgamma,gamma_var);
    
    currentlmux1 = calc_lmu(Za,currentgamma);
    
    currentlambda = sample_lambda(y1,currentx1,currentlambda,a0l,b0l,propl1,propl2);
    
    //std::cout << "6d\n";
    
    if((i>99)&&(i<burn)&&(i%20==0)){
      gamma_var = cov(gamma.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(gamma_var(j,j)==0){
          std::cout << "gamma proposal covariance matrix has variance 0\n";
        }
      }
    }
    if((i>99)&&(i<burn)&&(i%20==0)){
      betax_var = cov(betax.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(betax_var(j,j)==0){
          std::cout << "betax proposal covariance matrix has variance 0\n";
        }
      }
    }
    // if(i==burn){
    //   std::cout << "sigma betax proposal = " << betax_var << "\n";
    // }
    
    //std::cout << "7\n";
    
    currentsigma2y = sample_sigma2(a0y,b0y,y2sub1,y2sub2,x1x2p1,x1x2p2,currentbetay);
    //currentsigma2x = sample_sigma2(a0x,b0x,currentx2,Zbx,currentbetax);
    
    //std::cout << "8\n";
    
    // currentx1 = sample_x1(currentx1,currentx2,y1,y2,Zbx,currentp,exp(currentlmux1),
    //   currentlmux2,currentlmuy,currentbetay,currentbetax,
    //   currentalpha,currentsigma2x,currentsigma2y,currenteta,
    //   x1propa,x1propb,currentlambda);
    currentx1 = sample_x1(currentx1,currentx2,y1,y2,Zbx,currentp,exp(currentlmux1),
                          exp(currentlmux2),currentlmuy,currentbetay,currentbetax,
                          currentalpha,currenttheta,currentsigma2y,currenteta,
                          x1propa,x1propb,currentlambda,x1x2);
    
    x1x2 = arma::join_rows(intercept,arma::join_rows(currentx1,currentx2));
    x1x2y = arma::join_rows(intercept,arma::join_rows(currentx1,log(currentx2)));//x1x2;//arma::join_rows(intercept,currentx2);
    Zbx = arma::join_rows(Zb,currentx1);
    currentlmux2 = calc_lmu(Zbx,currentbetax);
    currentp = calc_p(x1x2,currentalpha);
    currentlmuy = calc_lmu(x1x2y,currentbetay);
    
    //std::cout << "9\n";
    
    // currentx2 = sample_x2(currentx1,currentx2,y2,currentp,currentlmux2,
    //                   currentlmuy,currentbetay,currentalpha,currentsigma2x,
    //                   currentsigma2y,propx2,vx2);
    currentx2 = sample_x2(currentx1,currentx2,y2,currentp,exp(currentlmux2),
                          currentlmuy,currentbetay,currentalpha,currenttheta,
                          currentsigma2y,propx2,vx2,x1x2);

    if((i>99)&&(i<burn)&&(i%20==0)){
      for(int j=0;j<n;j++){
        vx2[j] = var(latentx2.col(j).rows(0,i-1));
        if(vx2[j]==0){
          std::cout << "x2 proposal variance is 0 for individual " << j << " \n";
        }
      }
    }
    
    x1x2 = arma::join_rows(intercept,arma::join_rows(currentx1,currentx2));
    x1x2y = arma::join_rows(intercept,arma::join_rows(currentx1,log(currentx2)));//x1x2;//arma::join_rows(intercept,currentx2);
    
    x1x2p1 = subset2(x1x2y,indy1);
    x1x2p2 = subset2(x1x2y,indy2);
    x1x2p = arma::join_cols(x1x2p1,x1x2p2);
    
    //std::cout << "10\n";
    
    
    betay.row(i) = currentbetay.t();
    betax.row(i) = currentbetax.t();
    gamma.row(i) = currentgamma.t();
    alpha.row(i) = currentalpha.t();
    sigma2x[i] = currentsigma2x;
    sigma2y[i] = currentsigma2y;
    eta[i]      = currenteta;
    theta[i]      = currenttheta;
    lambda[i]      = currentlambda;
    delta[i]      = currentdelta;
    latentx1.row(i) = currentx1.t();
    latentx2.row(i) = currentx2.t();
    muy.row(i) = (currentlmuy.t());
    mux1.row(i) = exp(currentlmux1.t());
    mux2.row(i) = exp(currentlmux2.t());
    pg0.row(i) = currentp.t();
    
    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  } 
  return List::create(
    Named("betay") = betay,
    Named("betax") = betax,
    Named("alpha") = alpha,
    Named("gamma") = gamma,
    Named("sigma2x") = sigma2x,
    Named("sigma2y") = sigma2y,
    Named("eta")    = eta,
    Named("theta")    = theta,
    Named("delta")    = delta,
    Named("lambda")    = lambda,
    Named("latentx1") = latentx1,
    Named("latentx2") = latentx2,
    Named("muy") = muy,
    Named("mux1") = mux1,
    Named("mux2") = mux2,
    Named("p") = pg0);
} 








// !!!!!!!!!!!!!!

