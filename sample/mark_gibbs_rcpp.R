library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

one <- function() 1L
cppFunction('int one() {
  return 1;
}')
  
###### in R #######
MJgibbs_r <- function(nsamples = 5000){
  # Simulate data
  x1 <- rnorm(100,0,1)
  x2 <- rbinom(100,size=1,0.5)
  y <- 1 + x1 + 2*x2 + rnorm(100,0,1)
  
  n <- length(y)
  
  # Design matrix
  X <- as.matrix(cbind(rep(1,100),x1,x2))
  return(X)
}


one <- function() 1L
cppFunction('int one() {
            return 1;
            }')

cppFunction('arma::mat MJgibbs_cpp(int nsamples) {
  NumericVector x1 = rnorm(100,0,1);
  NumericVector x2 = rbinom(100,1,0.5);
  arma::vec y = 1 + x1 + 2*x2 + rnorm(100,0,1);
  int n = y.size();
  

  arma::mat X(100,3);
 int nrow = X.n_rows;
for(int i = 0; i < nrow; i++){
X(i,1) = x1[i];
X(i,2) = x2[i];
}
int p = X.n_cols;

double a  = 0.01;
double b  = 0.01;

double sigma_beta = 1e10;

arma::vec beta(p);
double sigma_epsilon = 1;

arma::mat post_samples(nsamples,4);

arma::mat sigma_BETA;
arma::mat mu_BETA;
arma::mat diag_n = arma::eye(n,n);
arma::mat diag_p = arma::eye(p,p);
arma::mat transX = arma::trans(X);
arma::mat b_star;
double a_star = a + n/2;
sigma_BETA = inv(transX * inv(sigma_epsilon * diag_n)* X + inv(sigma_beta * diag_p)); 

mu_BETA = sigma_BETA * transX * inv(sigma_epsilon * diag_n) * y;
arma::mat mvrnormArma();
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
arma::mat Y = arma::randn(n, ncols);
return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

for(int i =0; i < nsamples ; i++){
beta = mvrnormArma(1,mu_BETA,sigma_BETA);
b_star = trans(y - X * beta) * (y -X * beta)/2 + b;
}
return mu_BETA;      

}',depends="RcppArmadillo")
b_star <- t(y - X%*%beta)%*%(y - X%*%beta)/2 + b
sigma_EPSILON <- rigamma(1,a_star,b_star)
post.samples[i,] <- c(beta,sigma_EPSILON)
a_star <- a + n/2
  X <- as.matrix(cbind(rep(1,100),x1,x2))
  p <- dim(X)[2]
  
  # Prior specification
  # Hyperparameters for IG(a,b) prior for sigma_epsilon
  a <- b <- 0.01
  
  # prior variance for beta
  sigma_beta <- 1e10
  
  # Initial Values
  beta <- numeric(p)
  sigma_epsilon <- 1
  
  # Number of MCMC draws
  #nsamples <- 5000
  
  # Create matrix to store results
  post.samples <- matrix(nrow=nsamples,ncol=4)
}