#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
namespace rdist {
arma::mat mvrnormArma(int n, arma::mat mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

}

using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat MJgibbs_cpp(int nsamples) {
  NumericVector x1 = rnorm(100,0,1);
  NumericVector x2 = rbinom(100,1,0.5);
  arma::vec y = 1 + x1 + 2*x2 + rnorm(100,0,1);
  int n = y.size();
  
  
  arma::mat X(100,3);
  int nrow = X.n_rows;
  for(int i = 0; i < nrow; i++){
    X(i,0) = 1;
    X(i,1) = x1[i];
    X(i,2) = x2[i];
  }
  int p = X.n_cols;
  
  double a  = 0.01;
  double b  = 0.01;
  
  double sigma_beta = 1e10;
  
  arma::mat beta = arma::mat(p,1);
  /*double sigma_epsilon = 1;*/
  
  arma::mat post_samples(nsamples,4);
  
  arma::mat sigma_BETA;
  arma::mat mu_BETA;
  arma::mat diag_n = arma::eye(n,n);
  arma::mat diag_p = arma::eye(p,p);
  arma::mat transX = arma::trans(X);
  double a_star = a + n/2;
  sigma_BETA = inv(transX*X + inv(sigma_beta * diag_p)); 
  
  mu_BETA = sigma_BETA * transX* y;

  double b_star;
  /*double sigma_EPSILON;*/
  
  for(int i =0; i < nsamples ; i++){
    beta = rdist::mvrnormArma(1,mu_BETA,sigma_BETA);
    b_star= as_scalar(trans(y - (X * trans(beta))) * (y -(X * trans(beta)))/2 + b);
   /* sigma_EPSILON  = (1/(Rcpp::rgamma(1,a_star,b_star)));*/
    post_samples(i,0)=beta(0,0);
    post_samples(i,1)=beta(0,1);
    post_samples(i,2)=beta(0,2);
    post_samples(i,3) = 1/(R::rgamma(a_star,b_star));
 }
  return post_samples;
                             
  /* return mu_BETA;*/
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
system.time(MJgibbs_cpp(500))
*/
