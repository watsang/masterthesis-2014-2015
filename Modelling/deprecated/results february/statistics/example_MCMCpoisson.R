library(MCMCpack)
set.seed(11119)
n <- 150
x1 <- runif(n, 0, 0.5)
true.beta1 <- c(1,  1)
true.beta2 <- c(1,  -2)
true.beta3 <- c(1,  2)    

## set true two breaks at (50, 100)
true.s <- rep(1:3, each=n/3)
mu1 <- exp(1 + x1[true.s==1]*1)
mu2 <- exp(1 + x1[true.s==2]*-2)
mu3 <- exp(1 + x1[true.s==3]*2)

y <- as.ts(c(rpois(n/3, mu1), rpois(n/3, mu2), rpois(n/3, mu3)))
formula = y ~ x1

## fit multiple models with a varying number of breaks
model0 <-  MCMCpoissonChange(formula, m=0, 
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
model1 <-  MCMCpoissonChange(formula, m=1, 
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
model2 <-  MCMCpoissonChange(formula, m=2,
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
model3 <-  MCMCpoissonChange(formula, m=3, 
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
model4 <-  MCMCpoissonChange(formula, m=4, 
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    
model5 <-  MCMCpoissonChange(formula, m=5, 
                             mcmc = 1000, burnin = 1000, verbose = 500, 
                             b0 = rep(0, 2), B0 = 5*diag(2), marginal.likelihood = "Chib95")    

## find the most reasonable one
print(BayesFactor(model0, model1, model2, model3, model4, model5))

## draw plots using the "right" model
par(mfrow=c(attr(model2, "m") + 1, 1), mai=c(0.4, 0.6, 0.3, 0.05))
plotState(model2, legend.control = c(1, 0.6))
plotChangepoint(model2, verbose = TRUE, ylab="Density", start=1, overlay=TRUE)

## No covariate case
model2.1 <- MCMCpoissonChange(y ~ 1, m = 2, c0 = 2, d0 = 1,
                              mcmc = 1000, burnin = 1000, verbose = 500, 
                              marginal.likelihood = "Chib95")   
print(BayesFactor(model2, model2.1)) 
## End(Not run)

###########################################################################################################
###########################################################################################################
### N E G A T I V E    B I N O M I A L   R E G R E S S I O N
###########################################################################################################
###########################################################################################################

## Not run:
## logistic regression with an improper uniform prior
## X and y are passed as args to MCMCmetrop1R
logitfun <- function(beta, y, X){
  eta <- X %*% beta
  p <- 1.0/(1.0+exp(-eta))
  sum( y * log(p) + (1-y)*log(1-p) )
}
x1 <- rnorm(1000)
x2 <- rnorm(1000)
Xdata <- cbind(1,x1,x2)
p <- exp(.5 - x1 + x2)/(1+exp(.5 - x1 + x2))
yvector <- rbinom(1000, 1, p)

post.samp <- MCMCmetrop1R(logitfun, theta.init=c(0,0,0),
                          X=Xdata, y=yvector,
                          thin=1, mcmc=40000, burnin=500,
                          tune=c(1.5, 1.5, 1.5),
                          verbose=500, logfun=TRUE)
raftery.diag(post.samp)
plot(post.samp)
summary(post.samp)
## ##################################################
## negative binomial regression with an improper unform prior
## X and y are passed as args to MCMCmetrop1R
negbinfun <- function(theta, y, X){
  k <- length(theta)
  beta <- theta[1:(k-1)]
  alpha <- exp(theta[k])
  mu <- exp(X %*% beta)
  log.like <- sum(
    lgamma(y+alpha) - lfactorial(y) - lgamma(alpha) +
      alpha * log(alpha/(alpha+mu)) +
      y * log(mu/(alpha+mu))
  )
}
n <- 1000
x1 <- rnorm(n)
x2 <- rnorm(n)
XX <- cbind(1,x1,x2)
mu <- exp(1.5+x1+2*x2)*rgamma(n,1)
yy <- rpois(n, mu)
post.samp <- MCMCmetrop1R(negbinfun, theta.init=c(0,0,0,0), y=yy, X=XX,
                          thin=1, mcmc=35000, burnin=1000,
                          tune=1.5, verbose=500, logfun=TRUE,
                          seed=list(NA,1))
raftery.diag(post.samp)
plot(post.samp)
summary(post.samp)
## ##################################################
## sample from a univariate normal distribution with
## mean 5 and standard deviation 0.1
##
## (MCMC obviously not necessary here and this should
## really be done with the logdensity for better
## numerical accuracy-- this is just an illustration of how
## MCMCmetrop1R works with a density rather than logdensity)
post.samp <- MCMCmetrop1R(dnorm, theta.init=5.3, mean=5, sd=0.1,
                          thin=1, mcmc=50000, burnin=500,
                          tune=2.0, verbose=5000, logfun=FALSE)
summary(post.samp)