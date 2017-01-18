#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


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

// arma::vec calc_lmu(arma::mat Zb, arma::vec b, arma::vec beta){
//   int n = Zb.n_rows;
//   int k = Zb.n_cols;
//   arma::vec out = arma::zeros(n);
//   for(int i=0;i<n;i++){
//     for(int j=0;j<k;j++){
//       out[i] += Zb(i,j)*beta[j];
//     }
//     out[i] += b[i];
//   }
//   return out;
// }

// [[Rcpp::export]]
List ppred_ln(List sample, List data, int nreps) {
  std::cout << "a\n";
  arma::mat latentx         = as<arma::mat>(sample["latentx"]);
  std::cout << "a2\n";
  
  arma::mat beta            = as<arma::mat>(sample["beta"]);
  arma::mat alpha           = as<arma::mat>(sample["alpha"]);
  arma::mat b1              = as<arma::mat>(sample["b1"]);
  arma::mat b2              = as<arma::mat>(sample["b2"]);
  arma::vec gamma           = as<arma::vec>(sample["gamma"]);
  arma::vec sigma2          = as<arma::vec>(sample["sigma2"]);
  arma::vec sigma2b1        = as<arma::vec>(sample["sigma2b1"]);
  arma::vec sigma2b2        = as<arma::vec>(sample["sigma2b2"]);
  arma::vec sigma2w         = as<arma::vec>(sample["sigma2w"]);
  arma::vec sigma2y         = as<arma::vec>(sample["sigma2y"]);
  arma::vec corrb           = as<arma::vec>(sample["corrb"]);
  std::cout << "b\n";
  arma::mat Za              = as<arma::mat>(data["Za"]);
  arma::mat Zb              = as<arma::mat>(data["Zb"]);
  arma::mat y               = as<arma::mat>(data["y"]);
  arma::mat w               = as<arma::mat>(data["w"]);
  std::cout << "c\n";
  int n = Za.n_rows;
  std::cout << "d\n";
  //allocate storage
  arma::vec nzeros(nreps);
  std::cout << "e\n";
  //for calculations
  arma::vec pi(n);
  arma::vec lmu(n);
  arma::vec x;
  arma::mat b;
  double r;
  int zerocnt;
  
  for(int i=0;i<nreps;i++){
    std::cout << "1\n";
    b = arma::join_cols(b1.row(i),b2.row(i));
    std::cout << "2\n";
    pi = calc_p(Za,alpha.row(i).t(),b.t());
    std::cout << "3\n";
    //lmu = calc_lmu(Zb,b.row(1).t(),beta.row(i).t());
    std::cout << "4\n";
    zerocnt = 0;
    for(int j=0;j<n;j++){
      r = R::runif(0,1);
      if(r > pi[j]){
        x[j] = 0.0;
        zerocnt+=1;
      }
      else{
        //x[j] = R::rlnorm(lmu[j],pow(sigma2[i],0.5));
      }
    }
    
    nzeros[i] = zerocnt;
    
    //std::cout << i << "\n";
  }
  
  return List::create(
    Named("zeros") = nzeros);
}

