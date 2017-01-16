#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;




arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

const double log2pi = std::log(2.0 * M_PI);


// arma::vec dmvnrm_arma(arma::mat x,  
//                       arma::rowvec mean,  
//                       arma::mat sigma, 
//                       bool logd = false) { 
//   int n = x.n_rows;
//   int xdim = x.n_cols;
//   arma::vec out(n);
//   arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
//   double rootisum = arma::sum(log(rooti.diag()));
//   double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
//   
//   for (int i=0; i < n; i++) {
//     //std::cout << x.row(i) << "\n" << mean;
//     
//     arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
//     out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
//   }  
//   
//   if (logd == false) {
//     out = exp(out);
//   }
//   return(out);
// }

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

arma::vec subset(arma::vec x, arma::vec ind){
  int n = arma::accu(ind);
  arma::vec out(n);
  int count = 0;
  for(int i=0;i<n;i++){
    if(ind[i]==1){
      out[count] = x[i];
      count++;
    }
  }
  return out;
}

arma::mat subset(arma::mat x, arma::vec ind){
  int n = arma::accu(ind);
  int k = x.n_cols;
  arma::mat out(n,k);
  int count = 0;
  for(int i=0;i<n;i++){
    if(ind[i]==1){
      out.row(count) = x.row(i);
      count++;
    }
  }
  return out;
}

arma::vec sample_u(arma::mat Za, arma::vec x, arma::vec alpha, arma::mat b){
  int n = Za.n_rows;
  int k = Za.n_cols;
  arma::vec out(n);
  double check;
  double mu;
  for(int i=0;i<n;i++){
    mu = 0.0;
    for(int j=0;j<k;k++){
      mu += Za(i,j)*alpha[j];
    }
    check = R::rnorm(mu+b(i,0),1);
    if(x[i]==0){
      out[i] = -abs(check);
    }
    else{
      out[i] = abs(check);
    }
  }
  return out;
}

arma::vec sample_alpha(arma::mat Za, arma::mat b, arma::vec u, 
                       arma::vec mu0a, arma::mat V0a){
  int k = Za.n_cols;
  arma::vec out(k);
  arma::mat Va = arma::inv(arma::inv(V0a)+Za.t()*Za);
  arma::vec mua = Va*(arma::inv(V0a)*mu0a + Za.t()*(u-b.col(0)));
  
  out = (mvrnormArma(1,mua,Va).row(0).t());
  

  return out;
}

arma::vec calc_p(arma::mat Za, arma::vec alpha, arma::mat b){
  int n = Za.n_rows;
  int k = Za.n_cols;
  arma::vec out(n);
  double x;
  for(int i=0;i<n;i++){
    x=0.0;
    for(int j=0;j<k;j++){
      x += Za(i,j)*alpha(j);
    }
    out[i] = R::pnorm5(x+b(i,0),0,1,1,0);
  }
  return out;
}

arma::vec sample_beta(arma::mat Zbp, arma::vec b, arma::vec mu0b, arma::mat V0b,
                      double sigma2, arma::vec xp, arma::vec ind){
  int k = Zbp.n_cols;
  arma::vec out(k);
  arma::mat Vb = arma::inv(arma::inv(V0b)+Zbp.t()*Zbp/sigma2);
  arma::vec mub = Vb*(arma::inv(V0b)*mu0b + Zbp.t()*(arma::log(xp)-b)/sigma2);
  
  out = (mvrnormArma(1,mub,Vb).row(0).t());
  
  return out;
}

double sample_gamma(arma::vec x,double sigma2y,double V0g,double mu0g,arma::vec ybar){
  double out;
  int n = x.size();
  double total = 0.0;
  double total2 = 0.0;
  for(int i=0;i<n;i++){
    total+=pow(x[i],2);
    total2+=x[i]*ybar[i];
  }
  double Vg = 1/(total/sigma2y+1/V0g); 
  double mg = Vg*(mu0g/V0g+total2/sigma2y);
  out = R::rnorm(mg,Vg);
  
  return out;
}


arma::vec calc_lmu(arma::mat Zbp, arma::vec b, arma::vec beta){
  int n = Zbp.n_rows;
  int k = Zbp.n_cols;
  arma::vec out(n);
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      out[i] += Zbp(i,j)*beta[j];
    }
    out[i] += b[i];
  }
  return out;
}

double sample_sigma2(double a0, double b0, arma::vec xp, arma::mat Zbp,
                     arma::vec beta, arma::vec bp){
  double out;
  double a = xp.size()/2.0 + a0;
  double b = b0 + 0.5*arma::accu(pow(arma::log(xp)-Zbp*beta-bp,2.0));
  out = 1/R::rgamma(a,1/b);
  return(out);
}  

arma::mat sample_Sigmab(int nu0, arma::mat D0, arma::mat currentb){
  int n = currentb.n_rows;
  int df = nu0 + n;
  arma::mat scale = D0 + currentb.t()*currentb;
  arma::mat out(2,2);
  out = rinvwish(1,df,scale);
  return out;
}

double log_qb(double x, arma::vec b, double p, double mu, double sigma2, arma::mat Sigmab){
  double ll1 = dmvnrm_arma(b,arma::zeros(2),Sigmab,true);
  double ll2;
  if(x==0){
    ll2 = log(1-p);
  }
  else{
    ll2 = log(p) + R::dlnorm(x,mu,sqrt(sigma2),true);
  }
  return ll1+ll2;
}

double log_qx(double x, arma::rowvec w, arma::rowvec y, double p, double mu,
              double gamma, double sigma2, double sigma2w, double sigma2y){
  double ll1;
  double ll2=0.0;
  int k = w.size();
  
  if(x==0){
    ll1 = log(1-p);
  }
  else{
    ll1 = log(p) + R::dlnorm(x,mu,sqrt(sigma2),true);
  }
  for(int i=0;i<k;i++){
    ll2 += R::dnorm4(w[i],x,sqrt(sigma2w),true) + R::dnorm4(y[i],x*gamma,sqrt(sigma2y),true); 
  }
  return ll1+ll2;
}


double log_gx(double x, double p0){
  double ll;
  if(x==0){
    ll = log(1-p0);
  }
  else{
    ll = log(p0) + R::dexp(x,0.0095,true);
  }
  return ll;
}



arma::mat sample_b(arma::mat b, arma::cube rw_var, arma::mat Za, arma::mat Zb,
                   arma::vec alpha, arma::vec beta,arma::vec x,arma::vec p,
                   arma::vec mu, double sigma2, arma::mat Sigma2){
  int n = b.n_rows;
  int ka= Za.n_cols;
  int kb= Zb.n_cols;
  double xa;
  double xb;
  arma::vec propb(2);
  double propp;
  double propmu;
  double lacceptprob;
  arma::mat out(n,2);
  
  for(int i=0;i<n;i++){
    xa=0.0;
    xb=0.0;
    propb = mvrnormArma(1,b.row(i).t(),2.88*rw_var.slice(i)).row(0).t();
    for(int j=0;j<ka;j++){
      xa+=Za(i,j)*alpha[j];
    }
    for(int j=0;j<kb;j++){
      xb+=Zb(i,j)*beta[j];
    }
    propp = R::pnorm(xa+propb[0],0,1,1,0);
    propmu = xb + propb[1];
    
    lacceptprob = log_qb(x[i],propb,propp,propmu,sigma2,Sigma2)-
      log_qb(x[i],b.row(i).t(),p[i],mu[i],sigma2,Sigma2);
    if(log(R::runif(0,1))<lacceptprob){
      out.row(i) = propb;
    }
  }
  return out;
}
  
  
  /* calculates a block diagnoal covariance matrix, must have even number of cols!,
   * block diagnonals are 2x2 for EE and ES
   */
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

arma::mat update_tune(arma::mat b){
  int n = b.n_rows;
  arma::cube out(2,2,n);
  for(int i=0;i<n;i++){
    out.slice(i) = cpp_cov(b);
  }
  return out;
}


List sample_x(arma::vec x, arma::mat w, arma::mat y, arma::vec p, arma::vec mu,
              double gamma, double sigma2, double sigma2w, double sigma2y, double p0){
  int n = x.size();
  arma::vec out = x; 
  int prop0 = 0;
  int prop1 = 0;
  int accept0 = 0;
  int accept1 = 0;
  
  double r;
  double propx;
  double lacceptprob;
  
  for(int i=0;i<n;i++){
    r= R::runif(0,1);
    if(r<(1-p0)){
      propx = 0.0;
      prop0 += 1;
    }
    else{
      propx = R::rexp(0.0095);
      prop1 += 1;
    }
    lacceptprob = log_qx(propx,w.row(i),y.row(i),p[i],mu[i],gamma,sigma2,sigma2w,sigma2y) - 
        log_qx(x[i],w.row(i),y.row(i),p[i],mu[i],gamma,sigma2,sigma2w,sigma2y) -
        log_gx(propx,p0) + log_gx(x[i],p0);
    if(log(R::runif(0,1))<lacceptprob){
      out[i]=propx;
      if(propx==0){
        accept0 += 1;
      } 
      else{
        accept1 += 1;
      }
    }
  }
  return List::create(
    Named("x") = out,
    Named("prop0") = prop0,
    Named("prop1") = prop1,
    Named("accept0") = accept0,
    Named("accept1") = accept1);
}




List mcmc_me_ln(List data, 
                List init, 
                List priors, 
                const int nreps, 
                const int burn=1000,
                const double p0=0.85){
  
  //data
  arma::mat Za              = as<arma::mat>(data["Za"]);
  arma::mat Zb              = as<arma::mat>(data["Zb"]);
  arma::mat y               = as<arma::mat>(data["y"]);
  arma::mat w               = as<arma::mat>(data["w"]);
  
  int n                     = Za.n_rows;
  int nr                    = y.n_cols;
  int na                    = Za.n_cols;
  int nb                    = Zb.n_cols;
  arma::vec ybar            = mean(y,1);

  //starting values
  arma::vec currentbeta     = as<arma::vec>(init["currentbeta"]);
  arma::vec currentalpha    = as<arma::vec>(init["currentalpha"]);
  double currentsigma2      = as<double>(init["currentsigma2"]);
  double currentsigma2w     = as<double>(init["currentsigma2w"]);
  double currentsigma2y     = as<double>(init["currentsigma2y"]);
  double currentgamma       = as<double>(init["currentgamma"]);
  arma::mat currentb        = as<arma::mat>(init["currentb"]);
  arma::mat currentSigmab   = as<arma::mat>(init["currentSigmab"]);
  arma::vec currentx        = as<arma::vec>(init["currentx"]);
  arma::vec currentu        = arma::zeros(n);
  
  arma::vec tune            = as<arma::vec>(init["tune"]);
  arma::vec currentp(n);
  arma::vec currentlmu(n);
    
//priors
  arma::vec mu0a            = as<arma::vec>(priors["mu0a"]);
  arma::vec mu0b            = as<arma::vec>(priors["mu0b"]);
  double mu0g               = as<double>(priors["mu0g"]);
  arma::mat V0a             = as<arma::mat>(init["V0a"]);
  arma::mat V0b             = as<arma::mat>(init["V0b"]);
  double V0g                = as<double>(priors["V0g"]);
  double a0                 = as<double>(priors["a0"]);
  double b0                 = as<double>(priors["b0"]);
  double a0w                = as<double>(priors["a0w"]);
  double b0w                = as<double>(priors["b0w"]);
  double a0y                = as<double>(priors["a0y"]);
  double b0y                = as<double>(priors["b0y"]);
  int nu0                   = as<int>(priors["nu0"]);
  arma::mat D0              = as<arma::mat>(init["D0"]);
  
  
  //storage
  arma::mat beta(nreps,nb);
  arma::mat alpha(nreps,na);
  arma::vec gamma(nreps);
  arma::vec sigma2(nreps);
  arma::vec sigma2b1(nreps);
  arma::vec sigma2b2(nreps);
  arma::vec corrb(nreps);
  arma::vec sigma2w(nreps);
  arma::vec sigma2y(nreps);
  arma::mat b1(nreps,n);
  arma::mat b2(nreps,n);
  arma::mat latentx(nreps,n);
  
  //to check acceptance rates
  int prop0 = 0;
  int prop1 = 0;
  int accept0 = 0;
  int accept1 = 0;
  
  List holder;
  
  //
  arma::cube rw_var(2,2,n);
  for(int i=0;i<n;i++){
    rw_var.slice(i).diag() = tune;
    rw_var(0,1,i)=rw_var(1,0,i) = 0.0;
  }
  
  //values for functions
  arma::vec ind0 = check0(currentx);
  arma::mat currentxp = subset(currentx,ind0);
  arma::mat Zbp = subset(Zb,ind0);
  arma::mat currentbp = subset(currentb,ind0);
  
  double e2 = 0;
  double v2 = 0;

  
  
  for(int i=0;i<nreps;i++){
    
    // !!!!!!!!!!!!!!!!!!!!!!
    //where do i put  these? after sampling x probably
    ind0 = check0(currentx);
    currentxp = subset(currentx,ind0);
    Zbp = subset(Zb,ind0);
    currentbp = subset(currentb,ind0);
    // !!!!!!!!!!!!!!!!!!!!!!
    
    
    currentu = sample_u(Za, currentx, currentalpha, currentb);
    
    currentalpha = sample_alpha(Za,currentb,currentu,mu0a,V0a);
    currentp = calc_p(Za,currentalpha,currentb);
    
    currentbeta = sample_beta(Zbp,currentbp.col(1),mu0b,V0b,currentsigma2,currentxp,ind0);
    currentlmu = calc_lmu(Zbp,currentbp.col(1),currentbeta);
    
    currentsigma2 = sample_sigma2(a0,b0,currentxp,Zbp,currentbeta,currentbp.col(1));
    
    currentSigmab = sample_Sigmab(nu0,D0,currentb);
    
    
    currentb = sample_b(currentb,rw_var,Za,Zb,currentalpha,currentbeta,currentx,
                        currentp,currentlmu,currentsigma2,currentSigmab);
    
    if((i>20)&&(i<burn)&&(i%20==0)){
      for(int j=0;j<n;j++){
        rw_var.slice(j) = update_tune(arma::join_rows(b1.rows(0,i-1).col(j),b2.rows(0,i-1).col(j)));
      }
    }
    
    e2 = 0;
    v2 = 0;
    for(int ii=0;ii<nr;ii++){
      e2 += (accu(pow(y.col(ii)-currentgamma*currentx,2)));
      v2 += (accu(pow(w.col(ii)-currentx,2)));
    }
    currentsigma2y = 1/R::rgamma(a0y+nr*n/2.0,1.0/(b0y+e2/2.0));
    currentsigma2w = 1/R::rgamma(a0w+nr*n/2.0,1.0/(b0w+v2/2.0));

    
    currentgamma = sample_gamma(currentx,currentsigma2y,V0g,mu0g,ybar);
    
    holder = sample_x(currentx,w,y,currentp,currentlmu,currentgamma,
                      currentsigma2,currentsigma2w,currentsigma2y,p0);
    currentx = as<arma::vec>(holder["x"]);
    prop0 += as<int>(holder["prop0"]);
    prop1 += as<int>(holder["prop1"]);
    accept0 += as<int>(holder["accept0"]);
    accept1 += as<int>(holder["accept1"]);
    
    
    
    beta.row(i) = currentbeta.t();
    alpha.row(i) = currentalpha.t();
    gamma[i] = currentgamma;
    sigma2[i] = currentsigma2;
    sigma2b1[i] = currentSigmab(0,0);
    sigma2b2[i] = currentSigmab(1,1);
    corrb[i] = currentSigmab(0,1)/sqrt(currentSigmab(0,0)*currentSigmab(1,1));
    sigma2w[i] = currentsigma2w;
    sigma2y[i] = currentsigma2y;
    b1.row(i) = currentb.col(0).t();
    b2.row(i) = currentb.col(1).t();
    latentx.row(i) = currentx.t();
    
  } 
  return List::create(
    Named("beta") = beta,
    Named("alpha") = alpha,
    Named("gamma") = gamma,
    Named("sigma2") = sigma2,
    Named("sigma2b1") = sigma2b1,
    Named("sigma2b2") = sigma2b2,
    Named("sigma2w") = sigma2w,
    Named("sigma2y") = sigma2y,
    Named("corrb") = corrb,
    Named("b1") = b1,
    Named("b2") = b2,
    Named("latentx") = latentx,
    Named("prop0") = prop0,
    Named("prop1") = prop1,
    Named("accept0") = accept0,
    Named("accept1") = accept1);
  
  
} 