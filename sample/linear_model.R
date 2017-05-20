###############################################
# Simulate data and fit Bayesian Linear Model #
#

###############################################

# Housekeeping
#set.seed(467)
setwd("C:/Users/mark.janko/Dropbox/programming/BayesianLinearModel")

# Load Necessary Libraries
library(mvtnorm)
library(MASS)
library(pscl)
library(coda)

gibbs_m <- function(nsamples = 5000){
# Simulate data
x1 <- rnorm(100,0,1)
x2 <- rbinom(100,size=1,0.5)
y <- 1 + x1 + 2*x2 + rnorm(100,0,1)

n <- length(y)

# Design matrix
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
#nsamples <- nsamples

# Create matrix to store results
post.samples <- matrix(nrow=nsamples,ncol=4)

# Begin sampler
sigma_BETA <- solve(
  t(X)%*%solve(sigma_epsilon*diag(1,n))%*%X + solve(sigma_beta*diag(1,p)))
mu_BETA <- sigma_BETA%*%t(X)%*%solve(sigma_epsilon*diag(1,n))%*%y
a_star <- a + n/2
for (i in 1:nsamples){
    
    # Update Beta
        
        
        beta <- mvrnorm(1,mu_BETA,sigma_BETA)
                            
    # Update sigma_EPSILON
        
        b_star <- t(y - X%*%beta)%*%(y - X%*%beta)/2 + b
        sigma_EPSILON <- rigamma(1,a_star,b_star)
    
    # Store samples
    post.samples[i,] <- c(beta,sigma_EPSILON)
}
return(post.samples)
}
# End sampler

par(mfrow=c(2,2))
plot(post.samples[,1],type="l",ylab=expression(paste("posterior for ",beta[0])),xlab="iteration")
plot(post.samples[,2],type="l",ylab=expression(paste("posterior for ",beta[1])),xlab="iteration")
plot(post.samples[,3],type="l",ylab=expression(paste("posterior for ",beta[2])),xlab="iteration")
plot(post.samples[,4],type="l",ylab=expression(paste("posterior for ",sigma^2)),xlab="iteration")

# Compare to linear model via MLE
post.samples <- as.mcmc(post.samples)
summary(post.samples)

fit <- lm(y~x1 + x2)
coefs <- coef(fit) # Compare to mean or median (same for symmetric dist'ns like beta) of posterior samples
cis <- confint(fit) # coverage very similar to 2.5% and 97.5% quantiles of posterior distributions



save.image("LinearModel.RData")
