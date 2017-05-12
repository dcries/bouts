#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "rtn1.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// re in mean of y1 and in mean of y2, normal
// y2 weibull
//same as 6 but with power link instead of log link
// power of 2 and l
// no y1 in y2 regression, re correlated


//// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) { 
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma); 
}

// [[Rcpp::export]]
arma::vec rdirich(arma::vec a){
  int n = a.size();
  double total;
  arma::vec out(n);
  arma::vec gam(n);
  for(int i=0;i<n;i++){
    gam[i] = R::rgamma(a[i],1);
  }
  total = arma::accu(gam);
  out = gam/total;
  return out;
}

//// [[Rcpp::export]]
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
    ////std::cout << x.row(i) << "\n" << mean;
    
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

//// [[Rcpp::export]]
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
  ////std::cout << x.row(i) << "\n" << mean;
  
  arma::vec z = rooti * arma::trans( x.t() - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;
  //}
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double dgenpois(double x, double mu, double lambda, bool logd=true, double m = 1.0){
  double out;
  out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
  
  if(!logd){
    out=exp(out);
  }
  return out;
}

arma::vec dgenpois(double x, arma::vec mu, double lambda, bool logd=true){
  int n = mu.size();
  arma::vec out(n);
  out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
  if(!logd){
    out=exp(out);
  }
  return out;
}

// //// [[Rcpp::export]]
// arma::vec dgenpois(double x, arma::vec mu, double lambda, bool logd=true, double m = std::numeric_limits<double>::infinity()){
//   int n = mu.size();
//   arma::vec out(n);
//   out = log(mu*(1-lambda)) + (x-1)*log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
//   for(int i=0;i<n;i++){
//     if((lambda < 0) & (x > m)){
//       out[i] = -1*std::numeric_limits<double>::infinity();
//     }
//   }
//   if(!logd){
//     out=exp(out);
//   }
//   return out;
// }
// 
// //// [[Rcpp::export]]
// arma::vec dgenpois(arma::vec x, arma::vec mu, double lambda, bool logd=true, double m = std::numeric_limits<double>::infinity()){
//   int n = mu.size();
//   arma::vec out(n);
//   out = log(mu*(1-lambda)) + (x-1)%log(mu-lambda*(mu-x)) - (mu-lambda*(mu-x)) - lgamma(x+1);
//   for(int i=0;i<n;i++){
//     if((lambda < 0) & (x[i] > m)){
//       out[i] = -1*std::numeric_limits<double>::infinity();
//     }
//   }
//   if(!logd){
//     out=exp(out);
//   }
//   return out;
// }
// 
// //// [[Rcpp::export]]
// arma::vec dgenpois(arma::vec x, arma::vec mu, arma::vec lambda, bool logd=true, double m = std::numeric_limits<double>::infinity()){
//   int n = mu.size();
//   arma::vec out(n);
//   out = log(mu%(1-lambda)) + (x-1)%log(mu-lambda%(mu-x)) - (mu-lambda%(mu-x)) - lgamma(x+1);
//   for(int i=0;i<n;i++){
//     if((lambda[i] < 0) & (x[i] > m)){
//       out[i] = -1*std::numeric_limits<double>::infinity();
//     }
//   }
//   if(!logd){
//     out=exp(out);
//   }
//   return out;
// }

//does currentx > 0
//// [[Rcpp::export]]
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

//// [[Rcpp::export]]
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


arma::vec calc_p(arma::vec x1, double lambda){
  int n = x1.size();
  arma::vec out = arma::zeros(n);
  out = 1.0-dgenpois(0,x1,lambda,false);
  return out;
}

double calc_p(double x1, double lambda){
  double out = 0.0;
  out = 1.0-dgenpois(0,x1,lambda,false);
  return out;
}

//// [[Rcpp::export]]
arma::vec sample_beta_2(arma::mat Zbp, arma::vec mu0b, arma::mat V0b,
                        double sigma2, arma::vec xp, arma::vec bp){
  int k = Zbp.n_cols;
  arma::vec out(k);
  arma::mat Vb = arma::inv(arma::inv(V0b)+Zbp.t()*Zbp/sigma2);
  arma::vec mub = Vb*(arma::inv(V0b)*mu0b + Zbp.t()*(log(xp)-bp)/sigma2);
  ////std::cout << "Vb=" << sigma2 <<"\n";
  out = (mvrnormArma(1,mub,Vb).row(0).t());
  
  return out;
}

double log_qbeta(arma::vec cx, arma::mat vx, double phi, arma::mat mu1, 
                 arma::mat y2, arma::vec betay){
  double out = 0.0;
  int n = y2.n_rows;
  int k = y2.n_cols;
  
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(y2(i,j)>0){
        //out += R::dgamma(y2[i],delta,1/(delta/mu1[i]),true);
        out += R::dweibull(y2(i,j),phi,mu1(i,j),true);
        //std::cout << "dweibull= " << R::dweibull(y2(i,j),phi,mu1(i,j),true) << "\n";
        //std::cout << "mu(i,j)= " << mu1(i,j) << "\n";
      }
    }
  }
  out += dmvnrm_arma(betay,cx.t(),vx,true);
  return(out);
}

arma::vec sample_betay(arma::vec cx, arma::mat vx, double phi, arma::mat mu1, 
                       arma::mat y2, arma::mat Za1, arma::mat Za2, arma::vec betay, 
                       arma::mat Sigma_1, arma::vec currentb, double divisor=1.0){
  int k = cx.size();
  int n = y2.n_rows;
  double lacceptprob;
  arma::vec out = betay;
  arma::mat mu1prop(n,2);
  arma::mat proposalm;
  arma::vec proposal(k);
  // //std::cout << "gamma= " << gamma << "\n";
  // //std::cout << "Sigma= " << Sigma_1 << "\n";
  proposalm = mvrnormArma(1,betay,pow(2.4,2)*Sigma_1/(k*divisor));
  
  for(int i=0;i<k;i++){
    proposal[i] = proposalm(0,i);
  }
  
  mu1prop.col(0) = pow(Za1*proposal+currentb,2.0);
  mu1prop.col(1) = pow(Za2*proposal+currentb,2.0);
  
  lacceptprob = log_qbeta(cx,vx,phi,mu1prop,y2,proposal)-log_qbeta(cx,vx,phi,mu1,y2,betay);
  //std::cout << lacceptprob << "\n";
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return(out);
}

//// [[Rcpp::export]]
double log_qgamma(arma::vec cx, arma::mat vx, double lambda, arma::vec mu1, 
                  arma::mat y1, arma::mat y2, arma::mat muy2, arma::vec p, double phi, arma::vec gam){
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
        //out += log(p[i]) + R::dlnorm(y2(i,j),muy2(i,j),sqrt(sigma2y),true);
        out += log(p[i]) + R::dweibull(y2(i,j),phi,muy2(i,j),true);
        
      }
    }
    //out += R::dgamma(x1[i],eta,1/(eta/mu1[i]),true);
    
  }
  out += dmvnrm_arma(gam,cx.t(),vx,true);
  return(out);
}

//// [[Rcpp::export]]
arma::vec sample_gamma(arma::vec cx, arma::mat vx, double lambda, arma::vec mu1, 
                       arma::mat y1, arma::mat y2, arma::mat muy2, arma::mat Za, 
                       arma::vec gamma, arma::vec p, double phi,
                       arma::mat Sigma_1, arma::mat currentb,double divisor=1.0){
  int k = cx.size();
  int n = y1.n_rows;
  double lacceptprob;
  arma::vec out = gamma;
  arma::vec mu1prop(n);
  arma::vec propp(n);
  arma::mat proposalm;
  arma::vec proposal(k);
  // //std::cout << "before mvrnorm \n";
  // //std::cout << "gamma= " << gamma << "\n";
  // //std::cout << "Sigma= " << Sigma_1 << "\n";
  proposalm = mvrnormArma(1,gamma,pow(2.4,2)*Sigma_1/(k*divisor));
  
  ////std::cout << "after mvrnorm\n";
  for(int i=0;i<k;i++){
    proposal[i] = proposalm(0,i);
  }
  
  mu1prop = exp(Za*proposal+currentb.col(0));
  propp = calc_p(mu1prop,lambda);
  
  lacceptprob = log_qgamma(cx,vx,lambda,mu1prop,y1,y2,muy2,propp,phi,proposal)-
    log_qgamma(cx,vx,lambda,mu1,y1,y2,muy2,p,phi,gamma);
  ////std::cout << "after dvrnorm\n";
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return(out);
}



//// [[Rcpp::export]]
arma::vec calc_lmu(arma::mat Zb, arma::vec beta, arma::vec currentb){
  int n = Zb.n_rows;
  int k = Zb.n_cols;
  arma::vec out = arma::zeros(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out[i] += Zb(i,j)*beta[j];
    }
    out[i] += currentb[i];
  }
  return out;
}

//// [[Rcpp::export]]
double sample_sigma2(double a0, double b0, arma::vec xp, arma::mat Zbp,
                     arma::vec beta, arma::vec currentb){
  double out;
  double a = xp.size()/2.0 + a0;
  double b = b0 + 0.5*arma::accu(pow(log(xp)-Zbp*beta-currentb,2.0));
  out = 1/R::rgamma(a,1/b);
  return(out);
}  



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


double sample_sigma2(double a0, double b0, arma::vec u){
  double out;
  double a = u.size()/2.0 + a0;
  double b = b0 + 0.5*arma::accu(pow(u,2.0));
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

//// [[Rcpp::export]]
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


double log_ep(double ad, double bd, double phi, arma::mat y2, arma::mat x2){
  double out = (ad-1)*log(phi) - phi*bd;
  int n = y2.n_rows;
  int k = y2.n_cols;
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      if(y2(i,j)>0){
        //out += R::dgamma(y2(i,j),delta,1.0/(delta/x2[i]),true);
        out += R::dweibull(y2(i,j),phi,x2(i,j),true);
        
      }
    }
  }
  return(out);
}

double sample_phi(double ad, double bd, double phi, arma::mat y2, arma::mat x2, double propa, double propb){
  double proposal;
  double lacceptprob;
  double out = phi;
  
  proposal = R::rgamma(propa,1.0/propb);
  lacceptprob = log_ep(ad,bd,proposal,y2,x2) - log_ep(ad,bd,phi,y2,x2) -
    R::dgamma(proposal,propa,1.0/propb,true) + R::dgamma(phi,propa,1.0/propb,true);
  
  if(log(R::runif(0,1)) < lacceptprob){
    out = proposal;
  }
  return(out);
}



double log_ql(arma::mat y1, arma::mat y2, arma::vec x1, arma::mat muy2,
              arma::vec p,double phi,double lambda, double a0l, double b0l){
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
        //out += log(p[i]) + R::dlnorm(y2(i,j),muy2(i,j),sqrt(sigma2y),true);
        out += log(p[i]) + R::dweibull(y2(i,j),phi,muy2(i,j),true);
      }
    }
    
  }
  out += (a0l-1)*log(lambda) + (1-b0l)*log(1-lambda);
  return out;
}

//// [[Rcpp::export]]
double sample_lambda(arma::mat y1, arma::mat y2, arma::vec x1, arma::mat muy2,
                     arma::vec p, double phi, double lambda, double a0l, double b0l,
                     double propl1, double propl2){
  double out = lambda;
  double proposal = R::rbeta(propl1,propl2);
  double lacceptprob;
  int n = x1.size();
  arma::vec pprop(n);
  pprop = calc_p(x1,proposal);
  
  lacceptprob = log_ql(y1,y2,x1,muy2,pprop,phi,proposal,a0l,b0l) - log_ql(y1,y2,x1,muy2,p,phi,lambda,a0l,b0l) -
    R::dbeta(proposal,propl1,propl2,true) + R::dbeta(lambda,propl1,propl2,true);
  
  if(log(R::runif(0,1))<lacceptprob){
    out = proposal;
  }
  return out;
}

//// [[Rcpp::export]]
arma::mat sample_Sigmab(int d0, arma::mat D0, arma::mat currentb){
  int n = currentb.n_rows;
  int df = d0 + n;
  int k = currentb.n_cols;
  arma::mat scale = D0 + currentb.t()*currentb;
  arma::mat out(k,k);
  out = rinvwish(1,df,scale);
  return out;
}

double log_qb(arma::rowvec y1, arma::rowvec y2, double muy, arma::vec muy2, double lambdau,
              double phi, arma::vec currentb, arma::mat Sigmab){
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
      //ll += log(p0) + R::dlnorm(y2[i],muy2[i],sqrt(sigma2y),true);
      ll += log(p0) + R::dweibull(y2[i],phi,muy2[i],true);
      
    }
  }
  
  //ll += R::dgamma(x1,eta,1.0/(eta/mux1),true);
  //ll += R::dgamma(x2,theta,1.0/(theta/mux2),true);
  ////std::cout << dmvnrm_arma(currentb,zs,Sigmab,true) << "\n";
  ll += arma::accu(dmvnrm_arma(currentb,zs,Sigmab,true));
  //}
  
  return ll;
}




arma::mat sample_b(arma::mat y1, arma::mat y2,arma::vec muy1,arma::mat muy2,arma::mat currentb, 
                   arma::mat Sigmab, arma::vec gamma, arma::vec betay,
                   arma::mat Zg, arma::mat Zby1, arma::mat Zby2, double lambda, double phi, 
                   arma::cube Sigma_1){
  //std::cout << "a\n";
  
  int n = y1.n_rows;
  //int k = y2.n_cols;
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
    //std::cout << Sigma_1.slice(i) << "\n";
    propb2 = mvrnormArma(1,currentb.row(i).t(),pow(2.4,2)*Sigma_1.slice(i)/kp);
    //std::cout << "c\n";
    
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
    
    //std::cout << "d\n";
    
    lacceptprob = log_qb(y1.row(i),y2.row(i),exp(muy1prop),pow(muy2prop,2.0),lambda,phi,propb,Sigmab) - 
      log_qb(y1.row(i),y2.row(i),muy1[i],muy2.row(i).t(),lambda,phi,currentb.row(i).t(),Sigmab);
    //std::cout << "e\n";
    if(log(R::runif(0,1)) < lacceptprob){
      out.row(i) = propb.t();
    }
    
  }
  
  return out;
}

double log_qb1(arma::rowvec y1, arma::rowvec y2, double muy, arma::vec muy2, double lambdau,
               double phi, double currentb, double Sigmab){
  int n = y1.size();
  double ll=0.0;
  double p0;
  
  
  p0 = 1-dgenpois(0,muy,lambdau,false);
  for(int i=0;i<n;i++){
    ll += dgenpois(y1[i],muy,lambdau,true); 
    if(y2[i]==0){
      ll += log(1-p0);
    }
    else{
      //ll += log(p0) + R::dlnorm(y2[i],muy2[i],sqrt(sigma2y),true);
      ll += log(p0) + R::dweibull(y2[i],phi,muy2[i],true);
      
    }
  }
  
  //ll += R::dgamma(x1,eta,1.0/(eta/mux1),true);
  //ll += R::dgamma(x2,theta,1.0/(theta/mux2),true);
  ////std::cout << dmvnrm_arma(currentb,zs,Sigmab,true) << "\n";
  ll += R::dnorm(currentb,0,sqrt(Sigmab),true);
  //}
  
  return ll;
}

double log_qb2(arma::rowvec y2, arma::vec muy, double phi,
               double currentb, double Sigmab, double p0){
  int n = y2.size();
  double ll=0.0;
  
  for(int i=0;i<n;i++){
    if(y2[i]==0){
      ll += log(1-p0);
    }
    else{
      //ll += log(p0) + R::dlnorm(y2[i],muy2[i],sqrt(sigma2y),true);
      ll += log(p0) + R::dweibull(y2[i],phi,muy[i],true);
      
    }
    
  }
  
  ll += R::dnorm(currentb,0,sqrt(Sigmab),true); 
  //}
  
  return ll;
}




arma::vec sample_b1(arma::mat y1, arma::mat y2,arma::vec muy1,arma::mat muy2,arma::vec currentb, 
                    double Sigmab, arma::vec gamma, arma::vec betay,
                    arma::mat Zg, arma::mat Zby1, arma::mat Zby2, double lambda, double phi, 
                    arma::vec Sigma_1){
  //std::cout << "a\n";
  
  int n = y1.n_rows;
  //int k = y2.n_cols;
  int kp = currentb.n_cols;
  int kzg = Zg.n_cols;
  //int kzb = Zby1.n_cols;
  double lacceptprob;
  arma::vec out = currentb;
  double propb;
  double muy1prop;
  
  for(int i=0;i<n;i++){
    muy1prop = 0.0;
    //std::cout << "b\n";
    propb = R::rnorm(currentb[i],2.4*sqrt(Sigma_1[i]));
    //std::cout << "c\n";
    
    
    for(int j=0;j<kzg;j++){
      muy1prop += Zg(i,j)*gamma[j];
    }
    muy1prop += propb;
    
    //std::cout << "d\n";
    
    lacceptprob = log_qb1(y1.row(i),y2.row(i),exp(muy1prop),muy2.row(i).t(),lambda,phi,propb,Sigmab) - 
      log_qb1(y1.row(i),y2.row(i),muy1[i],muy2.row(i).t(),lambda,phi,currentb[i],Sigmab);
    //std::cout << "e\n";
    if(log(R::runif(0,1)) < lacceptprob){
      out[i] = propb;
    }
    
  }
  
  return out;
}

arma::vec sample_b2(arma::mat y2,arma::mat muy2,arma::vec currentb, 
                    double Sigmab, arma::vec betay,
                    arma::mat Zby1, arma::mat Zby2, double phi, 
                    arma::vec Sigma_1, arma::vec p0){
  //std::cout << "a\n";
  
  int n = y2.n_rows;
  //int k = y2.n_cols;
  int kzb = Zby1.n_cols;
  double lacceptprob;
  arma::vec out = currentb;
  double propb;
  arma::vec muy2prop(2);
  
  for(int i=0;i<n;i++){
    muy2prop = arma::zeros(2);
    //std::cout << "b\n";
    propb = R::rnorm(currentb[i],2.4*sqrt(Sigma_1[i]));
    //std::cout << "c\n";
    
    for(int j=0;j<kzb;j++){
      muy2prop[0] += Zby1(i,j)*betay[j];
      muy2prop[1] += Zby2(i,j)*betay[j];
      
    }
    muy2prop[0] += propb;
    muy2prop[1] += propb;
    
    //std::cout << "d\n";
    
    lacceptprob = log_qb2(y2.row(i),pow(muy2prop,2.0),phi,propb,Sigmab,p0[i]) - 
      log_qb2(y2.row(i),muy2.row(i).t(),phi,currentb[i],Sigmab,p0[i]);
    //std::cout << "e\n";
    if(log(R::runif(0,1)) < lacceptprob){
      out[i] = propb;
    }
    
  }
  
  return out;
}


// [[Rcpp::export]]
double calc_ind_dic(arma::mat y1, arma::mat y2, arma::vec mux1, arma::mat b,
                    double lambda, double sigma2y, arma::vec betay, 
                    arma::mat Sigma, arma::mat yZ1, arma::mat yZ2){
  int n = y1.n_rows;
  int k = y1.n_cols;
  int p = yZ1.n_cols;
  int nb = b.n_cols;
  double ll = 0.0;
  double ll2 = 0.0;
  double out;
  double prob;
  arma::vec muy1(k);
  arma::vec zs = arma::zeros(nb);
  arma::vec currb(nb);
  
  for(int i=0;i<n;i++){
    muy1 = arma::zeros(k);
    for(int p1=0;p1<p;p1++){
      muy1[0] += yZ1(i,p1)*betay[p1];
      muy1[1] += yZ2(i,p1)*betay[p1];
    }
    muy1[0] += b(i,1);
    muy1[1] += b(i,1);
    
    prob = (1-dgenpois(0,mux1[i],lambda,false));
    //std::cout << "i = " << i <<  " prob = " << prob << "\n";
    //prob=0.5;
    //std::cout << "e(muy1) " << exp(muy1) << "\n";
    for(int j=0;j<k;j++){
      ll2 += dgenpois(y1(i,j),mux1[i],lambda,true);
      //std::cout << "dgenpois = "  << dgenpois(y1(i,j),mux1[i],lambda+b(i,1),true,M[i]) << "\n";
      // std::cout << "ll after degenpois = " << ll << "\n";
      // std::cout << "mux1 = " << mux1[i] << "\n";
      // std::cout << "lambda + " << lambda+b(i,1) << "\n";
      if(y2(i,j)==0){
        ll2 += log(1-prob);
        //std::cout << "dlnorm = " << log(1-prob) << "\n";
        
      }
      else{
        ll2 += log(prob) + R::dlnorm(y2(i,j),muy1[j],sqrt(sigma2y),true);
        //std::cout << "dlnorm = " << log(prob) + R::dlnorm(y2(i,j),muy1[j],sqrt(sigma2y),true) << "\n";
        
      } 
    }
    //std::cout << "ll after dlnorm = " << ll << "\n";
    
    currb = b.row(i).t();
    ll += dmvnrm_arma(currb,zs.t(),Sigma,true);
    //std::cout << "dmvnrm_arma = "  << dmvnrm_arma(currb,zs.t(),Sigma,true) << "\n";
    
    //std::cout << "ll after dmvnorm = " << ll << "\n";
    
  }
  out = (ll2/k+ll) / n;
  return out;
}


// [[Rcpp::export]]
List mcmc_2part_nci7c(List data, 
                     List init, 
                     List priors, 
                     const int nreps, 
                     const int burn=1000){
  
  //const double infty = std::numeric_limits<double>::infinity();
  //std::numeric_limits<double>::infinity()
  
  //data
  arma::mat Za              = as<arma::mat>(data["Za"]);
  arma::mat Zb              = as<arma::mat>(data["Zb"]);
  arma::mat y1               = as<arma::mat>(data["y1"]);
  arma::mat y2               = as<arma::mat>(data["y2"]);
  
  //std::cout << "1\n";
  
  int n                     = Za.n_rows;
  //int nr                    = y1.n_cols;
  int na                    = Za.n_cols;
  int nb                    = Zb.n_cols;
  arma::vec ybar1            = mean(y1,1);
  arma::vec ybar2            = mean(y2,1);
  arma::vec intercept       = arma::ones(n);
  
  //std::cout << "2\n";
  
  
  //starting values
  arma::mat currentb = as<arma::mat>(init["currentb"]);
  arma::mat currentSigmab = as<arma::mat>(init["currentSigmab"]);
  arma::vec currentbetay     = as<arma::vec>(init["currentbetay"]);
  arma::vec currentgamma    = as<arma::vec>(init["currentgamma"]);
  arma::vec currentsigma2b(2);
  currentsigma2b[0] = currentsigma2b[1] = as<double>(init["currentsigma2b"]);
  
  double currentsigma2x     = as<double>(init["currentsigma2x"]);
  double currentsigma2y     = as<double>(init["currentsigma2y"]);
  double currenteta     = as<double>(init["currenteta"]);
  double currentlambda     = as<double>(init["currentlambda"]);
  double currentphi     = as<double>(init["currentphi"]);
  
  //std::cout << "2a\n";
  
  //double propa     = as<double>(init["propa"]);
  //double propb     = as<double>(init["propb"]);
  //double propax2     = as<double>(init["propax2"]);
  //double propbx2     = as<double>(init["propbx2"]);
  //double propx2     = as<double>(init["propx2"]);
  double propl1     = as<double>(init["propl1"]);
  double propl2     = as<double>(init["propl2"]);
  double propd1     = as<double>(init["propd1"]);
  double propd2     = as<double>(init["propd2"]);
  
  arma::vec currentx1        = as<arma::vec>(init["currentx1"]);
  arma::vec currentx2        = as<arma::vec>(init["currentx2"]);
  
  arma::vec gammatune            = as<arma::vec>(init["gammatune"]);
  arma::vec betaxtune            = as<arma::vec>(init["betaxtune"]);
  arma::vec btune            = as<arma::vec>(init["btune"]);
  
  arma::vec vx2             = as<arma::vec>(init["vx2"]);
  //arma::vec x1tune          = as<arma::vec>(init["x1tune"]);
  arma::vec x1propa          = as<arma::vec>(init["x1propa"]);
  arma::vec x1propb          = as<arma::vec>(init["x1propb"]);
  
  arma::vec currentp(n);
  arma::vec currentlmuy(n);
  arma::vec currentlmuy2(n);
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
  arma::mat D0             = as<arma::mat>(priors["D0"]);
  //double a0eta                 = as<double>(priors["a0eta"]);
  //double b0eta                 = as<double>(priors["b0eta"]);
  double a0delta                = as<double>(priors["a0delta"]);
  double b0delta                = as<double>(priors["b0delta"]);
  double a0x                = as<double>(priors["a0x"]);
  double b0x                = as<double>(priors["b0x"]);
  double a0y                = as<double>(priors["a0y"]);
  double b0y                = as<double>(priors["b0y"]);
  double a0l                = as<double>(priors["a0l"]);
  double b0l                = as<double>(priors["b0l"]);
  double d0                 = as<double>(priors["d0"]);
  //std::cout << "4\n";
  
  int nbre = currentb.n_cols;
  
  //storage
  arma::mat betay(nreps,mu0y2.size());
  arma::mat gamma(nreps,nb);
  arma::vec sigma2x(nreps);
  arma::vec sigma2y(nreps);
  arma::vec eta(nreps);
  arma::vec lambda(nreps);
  arma::mat latentx1(nreps,n);
  arma::mat latentx2(nreps,n);
  arma::mat muy(nreps,n);
  arma::mat mux1(nreps,n);
  arma::mat mux2(nreps,n);
  //arma::mat pg0(nreps,n);
  arma::mat b1(nreps,n);
  arma::mat b2(nreps,n);
  arma::mat sigma2b(nreps,nbre);
  arma::vec corrb(nreps);
  arma::mat m(nreps,n);
  arma::vec ind_dic(nreps-burn);
  arma::vec phi(nreps);
  double mean_dic;
  double dic;
  
  //std::cout << "5\n";
  
  arma::mat x1x2;
  arma::mat x1x2y1;
  arma::mat x1x2y2;
  arma::mat x1x2y;
  arma::mat x1x2p1;
  arma::mat x1x2p2;
  arma::mat x1x2p;
  
  arma::mat gamma_var = arma::zeros(na,na);
  gamma_var.diag() = gammatune;
  
  arma::mat betay_var = arma::zeros(nb,nb);
  betay_var.diag() = betaxtune;
  
  arma::cube b_var = arma::zeros(nbre,nbre,n);
  for(int i=0;i<n;i++){
    b_var.slice(i).diag() = btune;
  }
  // arma::mat b_var = arma::zeros(n,2);
  // for(int i=0;i<n;i++){
  //   b_var(i,0) = b_var(i,1) = btune[0];
  // }
  //int bvar0count = 0;
  //std::cout << "5\n";
  
  arma::vec indy1 = check0(y2.col(0));
  arma::vec indy2 = check0(y2.col(1));
  arma::vec indybar = check0(ybar2);
  
  arma::vec y2sub1 = subset(y2.col(0),indy1);
  arma::vec y2sub2 = subset(y2.col(1),indy2);
  arma::vec y2sub = arma::join_cols(y2sub1,y2sub2);
  
  //std::cout << "6\n";
  
  
  
  arma::mat Zbx;
  ////std::cout << "2\n";
  //std::cout << "7\n";
  Zbx = arma::join_rows(Zb,log(currentx1));
  x1x2 = Zb;//arma::join_rows(intercept,arma::join_rows(currentx1,currentx2));
  
  x1x2y1 = (Za);//arma::join_rows(Za,sqrt(y1.col(0)));//x1x2;//arma::join_rows(intercept,currentx2);
  x1x2y2 = (Za);//arma::join_rows(Za,sqrt(y1.col(1)));//x1x2;//arma::join_rows(intercept,currentx2);
  x1x2y = arma::join_rows(Za,ybar1);
  x1x2p1 = subset2(x1x2y1,indy1);
  x1x2p2 = subset2(x1x2y2,indy2);
  x1x2p = arma::join_cols(x1x2p1,x1x2p2);
  
  arma::mat b;
  arma::mat meanSigmab(nbre,nbre);
  
  arma::vec b3sub = arma::join_cols(subset(currentb.col(1),indy1),subset(currentb.col(1),indy2));
  currentlmuy = calc_lmu(x1x2y1,currentbetay,currentb.col(1));
  currentlmuy2 = calc_lmu(x1x2y2,currentbetay,currentb.col(1));
  
  int bvar0count;
  
  for(int i=0;i<nreps;i++){
    //std::cout << "4\n";
    
    //currenteta = sample_eta(a0eta,b0eta,currenteta,currentx1,exp(currentlmux1),propa,propb);
    //std::cout << "6a\n";
    //currentdelta = sample_delta(a0delta,b0delta,currentdelta,y2,currentx2,propd1,propd2);
    //currentbetay = sample_beta_2(x1x2p,mu0y2,V0y2,currentsigma2y,y2sub,b3sub);
    
    currentphi = sample_phi(a0delta,b0delta,currentphi,y2,arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),propd1,propd2);
    
    //std::cout << "5\n";
    
    currentbetay = sample_betay(mu0y2,V0y2,currentphi,arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
                                y2,x1x2y1,x1x2y2,currentbetay,betay_var,currentb.col(1),4.0);
    //std::cout << "5b\n";
    
    currentlmuy = calc_lmu(x1x2y1,currentbetay,currentb.col(1));
    currentlmuy2 = calc_lmu(x1x2y2,currentbetay,currentb.col(1));
    
    //std::cout << "6b\n";
    
    
    currentgamma = sample_gamma(mu0x1,V0x1,currentlambda,exp(currentlmux1),
                                y1,y2,arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
                                Za,currentgamma,currentp,currentphi,gamma_var,currentb,4.0);
    //std::cout << "6c\n";
    
    currentlmux1 = calc_lmu(Za,currentgamma,currentb.col(0));
    
    currentlambda = sample_lambda(y1,y2,exp(currentlmux1),arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
                                  currentp,currentphi,currentlambda,a0l,b0l,propl1,propl2);
    currentp = calc_p(exp(currentlmux1),currentlambda);
    //std::cout << "6d\n";
    
    //currentSigmab = sample_Sigmab(d0,D0,currentb);
    currentsigma2b[0] = sample_sigma2(a0x,b0x,currentb.col(0));
    currentsigma2b[1] = sample_sigma2(a0x,b0x,currentb.col(1));
    for(int j=0;j<2;j++){
      currentSigmab(j,j) = currentsigma2b[j];
    }
    //std::cout << "6d\n";
    
    if((i>99)&&(i<burn)&&(i%20==0)){
      betay_var = cov(betay.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(betay_var(j,j)==0){
          std::cout << "betax proposal covariance matrix has variance 0\n";
        }
      }
      //std::cout << "var = " << sqrt(betay_var.diag()) << "\n";
    }
    
    if((i>99)&&(i<burn)&&(i%20==0)){
      gamma_var = cov(gamma.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(gamma_var(j,j)==0){
          std::cout << "gamma proposal covariance matrix has variance 0\n";
        }
      }
    }
    
    currentb = sample_b(y1,y2,exp(currentlmux1),arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
                        currentb,currentSigmab,
                        currentgamma,currentbetay,Za,x1x2y1,x1x2y2,currentlambda,
                        currentphi,b_var);
    
    //std::cout << "1\n";
    // currentb.col(0) = sample_b1(y1,y2,exp(currentlmux1),arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
    //              currentb.col(0),currentsigma2b[0],
    //                                            currentgamma,currentbetay,Za,x1x2y1,x1x2y2,currentlambda,
    //                                            currentphi,b_var.col(0));
    //std::cout << "2\n";
    
    currentlmux1 = calc_lmu(Za,currentgamma,currentb.col(0));
    currentp = calc_p(exp(currentlmux1),currentlambda);
    //std::cout << "3\n";
    
    // currentb.col(1) = sample_b2(y2,arma::join_rows(pow(currentlmuy,2.0),pow(currentlmuy2,2.0)),
    //              currentb.col(1),currentsigma2b[1],currentbetay,x1x2y1,x1x2y2,currentphi,
    //              b_var.col(1),currentp);
    
    //std::cout << "7\n";
    
    // if((i>99)&&(i<burn)&&(i%20==0)){
    //   for(int j=0;j<n;j++){
    //     b = arma::join_rows(b1.col(j),b2.col(j));
    //     // b_var(j,0) = var(b1.col(j).rows(0,i-1));
    //     // b_var(j,1) = var(b2.col(j).rows(0,i-1));
    //     
    //     b_var.slice(j) = cov(b.rows(0,i-1));
    //   //   bvar0count = 0;
    //   //   if((b_var(0,0,j)==0) | (b_var(1,1,j)==0)){
    //   //          bvar0count += 1;
    //   //   }
    //   // }
    //   // if(bvar0count > 0){
    //   //       std::cout << "betax proposal covariance matrix has variance 0\n";
    //    }
    // }
    
    
    b3sub = arma::join_cols(subset(currentb.col(1),indy1),subset(currentb.col(1),indy2));
    currentx1 = exp(Za*currentgamma + currentb.col(0));
    //std::cout << "7b\n";
    
    //currentsigma2y = sample_sigma2(a0y,b0y,y2sub1,y2sub2,x1x2p1,x1x2p2,currentbetay,b3sub);
    
    //std::cout << "8\n";
    
    
    
    betay.row(i) = currentbetay.t();
    gamma.row(i) = currentgamma.t();    
    sigma2x[i] = currentsigma2x;
    sigma2y[i] = currentsigma2y;
    eta[i]      = currenteta;
    lambda[i]      = currentlambda;
    phi[i]      = currentphi;
    latentx1.row(i) = currentx1.t();
    latentx2.row(i) = currentx2.t();
    muy.row(i) = (currentlmuy.t());
    mux1.row(i) = exp(currentlmux1.t());
    mux2.row(i) = exp(currentlmux2.t());
    //pg0.row(i) = currentp.t();
    b1.row(i) = currentb.col(0).t();
    b2.row(i) = currentb.col(1).t();
    sigma2b.row(i) = (currentSigmab.diag()).t();
    corrb[i] = currentSigmab(0,1)/(sqrt(currentSigmab(0,0)*currentSigmab(1,1)));
    
    
    // if(i >= burn){
    //   ind_dic[i-burn] = calc_ind_dic(y1,y2,exp(currentlmux1),currentb,currentlambda,
    //                                  currentsigma2y,currentbetay,currentSigmab,x1x2y1,x1x2y2);
    //   //std::cout << "mean(ind_dic) = " << (ind_dic) << "\n";
    // }
    
    if(i % 1000==0){
      std::cout << "i= " << i << "\n";
    }
  } 
  
  // meanSigmab.diag() = (mean(sigma2b.rows(burn,nreps-1),0));
  // meanSigmab(0,1)=meanSigmab(1,0)=mean(corrb.subvec(burn,nreps-1))*sqrt(meanSigmab(0,0)*meanSigmab(1,1));
  // 
  // mean_dic = calc_ind_dic(y1,y2,trans(mean(mux1.rows(burn,nreps-1),0)),
  //                         arma::join_rows(trans(mean(b1.rows(burn,nreps-1),0)),trans(mean(b2.rows(burn,nreps-1),0))),
  //                         mean(lambda.subvec(burn,nreps-1)),mean(sigma2y.subvec(burn,nreps-1)),trans(mean(betay.rows(burn,nreps-1),0)),meanSigmab,x1x2y1,x1x2y2);
  // 
  // dic = -4*mean(ind_dic) + 2*mean_dic;
  //dic = 2*mean_dic;
  //std::cout << "mean(ind_dic) = " << (ind_dic) << "\n mean dic = " << mean_dic;
  
  return List::create(
    Named("betay") = betay,
    Named("gamma") = gamma,
    //Named("sigma2x") = sigma2x,
    Named("sigma2y") = sigma2y,
    Named("phi") = phi,
    //Named("eta")    = eta,
    Named("lambda")    = lambda,
    //Named("latentx1") = latentx1,
    //Named("latentx2") = latentx2,
    Named("b1") = b1,
    Named("b2") = b2,
    Named("sigma2b") = sigma2b,
    Named("corrb") = corrb,
    //Named("muy") = muy,
    //Named("mux1") = mux1,
    //Named("mux2") = mux2,
    Named("dic") = dic);
} 








// !!!!!!!!!!!!!!

