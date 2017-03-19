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


/*-------------------------------------------------------
# Generate Draws from an Inverse Wishart Distribution
# via the Bartlett Decomposition
#--------------------------------------------------------
# NOTE: output is identical to riwish from MCMCpack
#       provided the same random seed is used
#--------------------------------------------------------
#   n     number of samples
#   S     scale matrix 
#   v     degrees of freedom
#-------------------------------------------------------*/
// [[Rcpp::export]]
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
      x += Za(i,j)*alpha(j);
    }
    out[i] = R::pnorm5(x,0,1,1,0);
  }
  return out;
}

// [[Rcpp::export]]
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
                       arma::mat Sigma_1){
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
  proposalm = mvrnormArma(1,gamma,pow(2.4,2)*Sigma_1/k);
  
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

double sample_eta(double ae, double be, double currenteta, arma::vec x1, arma::vec mu1, double propa, double propb){
  double proposal;
  double lacceptprob;
  double out = currenteta;

  proposal = R::rgamma(propa,propb);
  lacceptprob = log_ex(ae,be,proposal,x1,mu1) - log_ex(ae,be,currenteta,x1,mu1) -
    R::dgamma(proposal,propa,propb,true) + R::dgamma(currenteta,propa,propb,true);
  
  if(log(R::runif(0,1)) < lacceptprob){
    out = proposal;
  }
  return(out);
}

//// [[Rcpp::export]]
// arma::mat sample_Sigmab(int nu0, arma::mat D0, arma::mat currentb){
//   int n = currentb.n_rows;
//   int df = nu0 + n;
//   arma::mat scale = D0 + currentb.t()*currentb;
//   arma::mat out(2,2);
//   out = rinvwish(1,df,scale);
//   return out;
// }

// // [[Rcpp::export]]
// double log_qb(double x, arma::rowvec w, arma::rowvec y, arma::vec b, double p, double mu, double sigma2,double sigma2w,double sigma2y, arma::mat Sigmab, double gamma){
//   arma::vec zs = arma::zeros(2);
//   int k = w.size();
//   
//   //arma::mat f(1,2);
//   //f.row(0) = arma::zeros(2);
//   //std::cout << "b= " << b << "\n zs = " << zs << "\n Sigmab = " << "\n";
//   double ll1 = dmvnrm_arma(b,zs.t(),Sigmab,true);
//   double ll2 = 0.0;
//   double ll3 = 0.0;
//   double ll4 = R::dlnorm(x,mu,sqrt(sigma2),true);
//   
//   for(int i=0;i<k;i++){
//     if(w[i]==0){
//       ll2 += log(1-p);
//     }
//     else{
//       ll2 += log(p) + R::dlnorm(w[i],log(x)/p,sqrt(sigma2w),true);
//     }
//     if(y[i]==0){
//       ll3 += log(1-p);
//     }
//     else{
//       ll3 += log(p) + R::dlnorm(y[i],gamma*log(x)/p,sqrt(sigma2y),true);
//     }
//   }
//   
//   
//   return ll1+ll2+ll3+ll4;
// }

// // [[Rcpp::export]]
double log_qx1(double x1, double x2, arma::rowvec y1, arma::rowvec y2, double p, double muy2,
              double mux1, double mux2,  double sigma2x, double sigma2y, double eta){
  double ll1;
  double ll2;
  double ll3=0.0;
  double ll4=0.0;
  
  int k = y2.size();
  
  ll1 = R::dlnorm(x2,mux2,sqrt(sigma2x),true);
  ll2 = R::dgamma(x1,eta,1/(eta/mux1),true);

    for(int i=0;i<k;i++){
      ll4 += R::dpois(y1[i],x1,true);
    if(y2[i]==0){
      ll3 += log(1-p);
    }
    else{
      ll3 += log(p) + R::dlnorm(y2[i],muy2,sqrt(sigma2y),true);
    }
  }
  
  return ll1+ll2+ll3+ll4;
}

double log_qx2(double x, arma::rowvec y, double p, double muy2,
               double mux2,  double sigma2x, double sigma2y){
  double ll1;
  double ll3=0.0;
  
  int k = y.size();
  
  ll1 = R::dlnorm(x,mux2,sqrt(sigma2x),true);
  
  
  for(int i=0;i<k;i++){
    if(y[i]==0){
      ll3 += log(1-p);
    }
    else{
      ll3 += log(p) + R::dlnorm(y[i],muy2,sqrt(sigma2y),true);
    }
  }
  
  return ll1+ll3;
}




/* calculates a block diagnoal covariance matrix, must have even number of cols!,
* block diagnonals are 2x2 for EE and ES
*/
// [[Rcpp::export]]
arma::mat cpp_cov(arma::mat x){
  int n = x.n_rows;
  int m = x.n_cols;
  arma::mat out = arma::zeros(m,m);
  arma::rowvec xbar = mean(x,0);
  double covar = 0.0;
  for(int i=0;i<m;i++){
    out(i,i) = var(x.col(i));
    if(i % 2 == 0){
      covar = 0.0;
      for(int j=0;j<n;j++){
        covar += (x(j,i)-xbar[i])*(x(j,i+1)-xbar[i+1]);
      }
      out(i,i+1) = out(i+1,i) = covar/(double(n)-1.0);
    }
  }
  return(out);
}



// // [[Rcpp::export]]
// arma::mat update_tune(arma::mat b){
//   int k = b.n_cols;
//   arma::mat out(k,k);
//   out = cpp_cov(b);
//   return out;
// }

// [[Rcpp::export]]
arma::vec sample_x2(arma::vec x1, arma::vec x2, arma::mat y2, arma::vec p, arma::vec mux2, 
                    arma::vec muy2, arma::vec betay, arma::vec alpha, double sigma2x, double sigma2y,
                    double propx2, arma::vec vx2){
  int n = x2.size();
  //int k = betay.size();
  double muyprop;
  double pprop;
  arma::vec out = x2; 
  
  
  double propx;
  double lacceptprob;
  
  for(int i=0;i<n;i++){
    //propx = R::rexp(propx2);
    propx = R::rnorm(x2[i],2.4*sqrt(vx2[i]));
    
    if(propx < 0.0){
      propx = abs(propx);
    }
    // muyprop = 0.0;
    // for(int j=0;j<k;j++){
    //   muyprop = 
    // }
    muyprop = betay[0] + betay[1]*x1[i] + betay[2]*propx;
    pprop = R::pnorm5(alpha[0] + alpha[1]*x1[i] + alpha[2]*propx,0,1,1,0);

    // lacceptprob = log_qx2(propx,y2.row(i),pprop,muyprop,mux2[i],sigma2x,sigma2y) -
    //   log_qx2(x2[i],y2.row(i),p[i],muy2[i],mux2[i],sigma2x,sigma2y) -
    //   R::dexp(propx,propx2,true) + R::dexp(x2[i],propx2,true);
    
    lacceptprob = log_qx2(propx,y2.row(i),pprop,muyprop,mux2[i],sigma2x,sigma2y) -
     log_qx2(x2[i],y2.row(i),p[i],muy2[i],mux2[i],sigma2x,sigma2y);
      
    if(log(R::runif(0,1))<lacceptprob){
      out[i]=propx;
    }
  }
  return out;
}

arma::vec sample_x1(arma::vec x1, arma::vec x2, arma::mat y1, arma::mat y2, 
                    arma::mat z, arma::vec p, arma::vec mux1, arma::vec mux2, 
                    arma::vec muy, arma::vec betay, arma::vec betax, 
                    arma::vec alpha, double sigma2x, double sigma2y,double eta){
  int n = x2.size();
  int k = z.n_cols;
  double muyprop;
  double pprop;
  double muxprop;
  double zb;
  arma::vec out = x1; 
  

  double propx;
  double lacceptprob;
  
  for(int i=0;i<n;i++){
    propx = R::rgamma(1,1);

    // muyprop = 0.0;
    // for(int j=0;j<k;j++){
    //   muyprop = 
    // }
    zb = 0.0;
    for(int j=2;j<k;j++){
      zb += z(i,j)*betax[j+2];
    }
    

    muyprop = betay[0] + betay[1]*x1[i] + betay[2]*propx;
    pprop = R::pnorm5(alpha[0] + alpha[1]*x1[i] + alpha[2]*propx,0,1,1,0);
    muxprop = betax[0] + betax[1]*x1[i] + zb;

    lacceptprob = log_qx1(propx,x2[i],y1.row(i),y2.row(i),pprop,muyprop,mux1[i],muxprop,sigma2x,sigma2y,eta) - 
      log_qx1(x1[i],x2[i],y1.row(i),y2.row(i),p[i],muy[i],mux1[i],mux2[i],sigma2x,sigma2y,eta) -
      R::dgamma(propx,1,1,true) + R::dgamma(x1[i],1,1,true);

    if(log(R::runif(0,1))<lacceptprob){
      out[i]=propx;
    }
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
  double propa     = as<double>(init["propa"]);
  double propb     = as<double>(init["propb"]);
  double propx2     = as<double>(init["propx2"]);
  
  arma::vec currentx1        = as<arma::vec>(init["currentx1"]);
  arma::vec currentx2        = as<arma::vec>(init["currentx2"]);
  arma::vec currentu        = arma::zeros(n);
  
  arma::vec gammatune            = as<arma::vec>(init["gammatune"]);
  arma::vec vx2             = as<arma::vec>(init["vx2"]);
  
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
  double a0x                = as<double>(priors["a0x"]);
  double b0x                = as<double>(priors["b0x"]);
  double a0y                = as<double>(priors["a0y"]);
  double b0y                = as<double>(priors["b0y"]);

  //std::cout << "4\n";
  
  //storage
  arma::mat betay(nreps,3);
  arma::mat betax(nreps,nb+1);
  arma::mat gamma(nreps,nb);
  arma::mat alpha(nreps,3);
  arma::vec sigma2x(nreps);
  arma::vec sigma2y(nreps);
  arma::vec eta(nreps);
  arma::mat latentx1(nreps,n);
  arma::mat latentx2(nreps,n);
  
  //std::cout << "5\n";
  
  arma::mat x1x2;
  arma::mat x1x2p1;
  arma::mat x1x2p2;
  arma::mat x1x2p;
  
  //std::cout << "1\n";
  //
  // arma::cube rw_var(2,2,n);
  // for(int i=0;i<n;i++){
  //   rw_var.slice(i).diag() = tune;
  //   rw_var(0,1,i)=rw_var(1,0,i) = 0.0;
  // }
  arma::mat gamma_var = arma::zeros(na,na);
  gamma_var.diag() = gammatune;
  
  //std::cout << "5\n";
  
  //values for functions
  // arma::vec ind0 = check0(currentx);
  // arma::vec currentxp = subset(currentx,ind0);
  // arma::mat Zbp = subset(Zb,ind0);
  // arma::mat currentbp = subset(currentb,ind0);

  arma::vec indy1 = check0(y2.col(0));
  arma::vec indy2 = check0(y2.col(1));
  arma::vec indybar = check0(ybar2);
  
  arma::vec y2sub1 = subset(y2.col(0),indy1);
  arma::vec y2sub2 = subset(y2.col(1),indy2);
  arma::vec y2sub = arma::join_cols(y2sub1,y2sub2);
  arma::vec ybarp2 = subset(ybar2,indybar);
  
  //std::cout << "6\n";
  
  arma::vec currentxy1 = subset(currentx1,indy1);
  arma::vec currentxy2 = subset(currentx1,indy2);
  arma::vec currentxybar = subset(currentx1,indybar);
  
  arma::vec currentpy1 = subset(currentp,indy1);
  arma::vec currentpy2 = subset(currentp,indy2);
  arma::vec currentpybar = subset(currentp,indybar);
  
  
  arma::mat Zbx;
  //std::cout << "2\n";
  //std::cout << "7\n";
  Zbx = arma::join_rows(Zb,currentx1);
  
  
  
  for(int i=0;i<nreps;i++){
    //std::cout << "4\n";
    
    x1x2 = arma::join_rows(intercept,arma::join_rows(currentx1,currentx2));
    //std::cout << "4a\n";
    
    x1x2p1 = subset2(x1x2,indy1);
    x1x2p2 = subset2(x1x2,indy2);
    //std::cout << "4b\n";
    
    x1x2p = arma::join_cols(x1x2p1,x1x2p2);
    //std::cout << "4c\n";
    
    // arma::mat Zbp = subset(Zb,ind0);
    
    currentu = sample_u(x1x2, ybar2, currentalpha);
    //std::cout << "5\n";
    
    currentalpha = sample_alpha(x1x2,currentu,mu0a,V0a);
    currentp = calc_p(x1x2,currentalpha);
    //std::cout << "6\n";
    
    currenteta = sample_eta(a0eta,b0eta,currenteta,currentx1,exp(currentlmux1),propa,propb);
    //std::cout << "6a\n";
    
    currentbetay = sample_beta_2(x1x2p,mu0y2,V0y2,currentsigma2y,y2sub);
    currentlmuy = calc_lmu(x1x2,currentbetay);
    
    //std::cout << "6b\n";
    
    //currentbetax = sample_beta_2(Zbx,mu0x2,V0x2,currentsigma2x,currentx2);
    currentlmux1 = calc_lmu(Zb,currentbetax);
    
    //std::cout << "6c\n";
    
    currentgamma = sample_gamma(mu0x1,V0x1,currenteta,exp(currentlmux1),currentx1,
                                Za,currentgamma,gamma_var);
    
    currentlmux2 = calc_lmu(Zbx,currentgamma);
    
    //std::cout << "6d\n";
    
    if((i>99)&&(i<burn)&&(i%20==0)){
      gamma_var = cov(gamma.rows(0,i-1));
      for(int j=0;j<na;j++){
        if(gamma_var(j,j)==0){
          std::cout << "gamma proposal covariance matrix has variance 0\n";
        }
      }
      //std::cout  <<"calc cov \n";
      
    }
    //std::cout << "6e\n";
    
    
    
    //std::cout << "7\n";
    
    currentsigma2y = sample_sigma2(a0y,b0y,y2sub1,y2sub2,x1x2p1,x1x2p2,currentbetay);
    currentsigma2x = sample_sigma2(a0x,b0x,currentx2,Zbx,currentbetax);
    
    //std::cout << "8\n";
    
    currentx1 = sample_x1(currentx1,currentx2,y1,y2,Zb,currentp,exp(currentlmux1),
                          currentlmux2,currentlmuy,currentbetay,currentbetax,
                          currentalpha,currentsigma2x,currentsigma2y,currenteta);
    
    Zbx = arma::join_rows(Zb,currentx1);
    currentlmux2 = calc_lmu(Zbx,currentgamma);
    
    //std::cout << "9\n";
    
    currentx2 = sample_x2(currentx1,currentx2,y2,currentp,currentlmux2,
                          currentlmuy,currentbetay,currentalpha,currentsigma2x,
                          currentsigma2y,propx2,vx2);

    if((i>99)&&(i<burn)&&(i%20==0)){
      for(int j=0;j<n;j++){
        vx2[j] = var(latentx2.col(j).rows(0,i-1));
        if(vx2[j]==0){
          std::cout << "x2 proposal variance is 0 for individual" << j << " \n";
        }
        
      }
      //std::cout  <<"calc cov \n";
    }
    //std::cout << "10\n";
    
    
    betay.row(i) = currentbetay.t();
    betax.row(i) = currentbetax.t();
    gamma.row(i) = currentgamma.t();
    alpha.row(i) = currentalpha.t();
    sigma2x[i] = currentsigma2x;
    sigma2y[i] = currentsigma2y;
    eta[i]      = currenteta;
    latentx1.row(i) = currentx1.t();
    latentx2.row(i) = currentx2.t();
    
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
    Named("latentx1") = latentx1,
    Named("latentx2") = latentx2);
} 








// !!!!!!!!!!!!!!

