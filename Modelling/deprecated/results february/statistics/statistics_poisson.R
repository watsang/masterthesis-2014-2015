library("rjags")

# Exercise births
# ---------------

setwd("D:/univ/2014-2015/thesis/KERMIT/results february")
bacteria <- read.csv("design_results_febr_counts.csv")

pathogen <- bacteria$LMG7866
evenness <- bacteria$evenness
cellcount <- bacteria$cellcount

# Simulate IQ data
set.seed(66)

n <- 20
x <- runif(n,-4,4)
epsilon <- rnorm(n,0,sigma)
y <- beta[1] + beta[2]*x + beta[3]*x^2 + epsilon

model <- lm(y~1+x+I(x^2))
plot(x,y)
sorted <- order(x)
model <- glm(pathogen~evenness+cellcount, family = quasipoisson)



# Frequentist analysis
summary(lm(births~1))

cat("\n      model\n      {\n        # priors\n        beta0 ~ dnorm(0,0.001)\n        beta1 ~ dnorm(0,0.001)\n        beta2 ~ dnorm(0,0.001)\n        \n        # likelihood\n        for(i in 1:N.cells)\n        {\n          n50[i] ~ dpois(lambda[i])\n          log(lambda[i]) <- beta0 + beta1*elev50[i] + beta2*pow(elev50[i],2)\n          # this part is here in order to make nice prediction curves:\n          prediction[i] ~ dpois(lambda[i])\n        } \n      }\n  ", 
    file = "model.txt")

modelstring <- "
model
{
  # priors
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  
  # likelihood
  for(i in 1:N.cells)
  {
    n50[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*elev50[i] + beta2*pow(elev50[i],2)
    # this part is here in order to make nice prediction curves:
    prediction[i] ~ dpois(lambda[i])
  } 
}"

# Bayesian analysis with JAGS
modelstring <- "
model{
  # Likelihood
  for (i in 1:nobs){
        y[i] ~ dnorm( mu[i] , tau )
        y[i] <- beta0 + beta1 * x[i] + beta2 * pow(x[i],2) 
    }
    beta0 ~ dnorm( 0 , 1.0E-12 )
    beta1 ~ dnorm( 0 , 1.0E-12 )
    beta2 ~ dnorm( 0 , 1.0E-12 ) 
    tau ~ dgamma( 0.001 , 0.001 )
#     y[i] ~ dnorm( mu[i] , tau )
#     mu[i] <- beta0 + beta1 * x[i] + beta2 * pow(x[i],2) # CHANGED
#     y[i] ~ beta0 + beta1*x[i] + beta2*pow(x[i],2)
#   }
  
  # Priors
#   beta0~dunif(-20,20)
#   beta1~dunif(-20,20)
#   beta2~dunif(-20,20)
  # level ~ dunif(0,1)
}
"

# Initiate the model
model <- jags.model(textConnection(modelstring),
                    data=list('y'=y,'nobs'=n),
                    #inits=list(births.boy=.5),
                    inits=list(list(beta0=0, beta1=0,beta2=0)),#,
                              #list(beta0=0, beta1=0,beta2=0),
                              #list(beta0=0, beta1=0,beta2=0)),
                    n.chains = 1,
                    n.adapt = 1000)


# Update initial model 
update(model, 1000)

# Run MCMC chain
out <- coda.samples(model,
                    variable.names=c('births.boy'),
                    n.iter=10000)

# Diagnostics checks
xyplot(out) # traceplots
gelman.plot(out) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0

# Results
plot(out)
summary(out)
HPDinterval(out)

# # --------------------------------------------------------------
# modelstring = "
# model {
#     for( i in 1 : Ndata ) {
# y[i] ~ dnorm( mu[i] , tau )
# mu[i] <- beta0 + beta1 * x[i] + beta2 * pow(x[i],2) # CHANGED
# }
# beta0 ~ dnorm( 0 , 1.0E-12 )
# beta1 ~ dnorm( 0 , 1.0E-12 )
# beta2 ~ dnorm( 0 , 1.0E-12 ) # CHANGED
# tau ~ dgamma( 0.001 , 0.001 )
# }
# " # close quote for modelstring
# writeLines(modelstring,con="model.txt")
# 
# #------------------------------------------------------------------------------
# # THE DATA.
# 
# # Generate some random data with a quadratic trend:
# set.seed(47405)
# N=50
# x = rnorm( N , mean=0 , sd=3 ) 
# y = rnorm( N )*0.4 + ( -4 + 2.6*(x-mean(x)) + 0.23*(x-mean(x))^2 )
# x=x+4
# y=y+15
# Ndata = length(y)
# 
# # Specify data, as a list.
# dataList = list(
#   x = x ,
#   y = y ,
#   Ndata = Ndata
# )
# 
# #------------------------------------------------------------------------------
# # INTIALIZE THE CHAINS.
# 
# # Use R's built-in least-squares regression to get plausible initial values:
# lmInfo = lm( dataList$y ~ dataList$x ) 
# b0Init = lmInfo$coef[1]
# bInit = lmInfo$coef[2]
# tauInit = length(dataList$y) / sum(lmInfo$res^2)
# initsList = list(
#   beta0 = b0Init ,
#   beta1 = bInit ,
#   beta2 = 0 , # CHANGED
#   tau = tauInit
# )
# 
# #------------------------------------------------------------------------------
# # RUN THE CHAINS
# 
# parameters = c("beta0" , "beta1" , "beta2" , "tau")  # CHANGED
# adaptSteps = 1000             # Number of steps to "tune" the samplers.
# burnInSteps = 2000            # Number of steps to "burn-in" the samplers.
# nChains = 3                   # Number of chains to run.
# numSavedSteps=100000           # Total number of steps in chains to save.
# thinSteps=1                   # Number of steps to "thin" (1=keep every step).
# nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# # Create, initialize, and adapt the model:
# jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
#                         n.chains=nChains , n.adapt=adaptSteps )
# # Burn-in:
# cat( "Burning in the MCMC chain...\n" )
# update( jagsModel , n.iter=burnInSteps )
# # The saved MCMC chain:
# cat( "Sampling final MCMC chain...\n" )
# codaSamples = coda.samples( jagsModel , variable.names=parameters , 
#                             n.iter=nPerChain , thin=thinSteps )
# # resulting codaSamples object has these indices: 
# #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
# 
# #------------------------------------------------------------------------------
# # EXAMINE THE RESULTS
# 
# checkConvergence = FALSE
# if ( checkConvergence ) {
#   openGraph(width=7,height=7)
#   autocorr.plot( codaSamples[[1]] , ask=FALSE ) 
#   show( gelman.diag( codaSamples ) )
#   effectiveChainLength = effectiveSize( codaSamples ) 
#   show( effectiveChainLength )
# }
# 
# # Convert coda-object codaSamples to matrix object for easier handling.
# # But note that this concatenates the different chains into one long chain.
# # Result is mcmcChain[ stepIdx , paramIdx ]
# mcmcChain = as.matrix( codaSamples )
# chainLength = NROW(mcmcChain)
# 
# # For convenience later, append a column with tau converted to sigma:
# sigma = 1 / sqrt( mcmcChain[, "tau" ] ) # Convert precision to SD
# mcmcChain = cbind( mcmcChain , sigma )
# 
# # Specify preferred file type for saved graphs:
# graphFileType = "jpg"
# 
# # Posterior prediction:
# # Specify x values for which predicted y's are wanted:
# extrapolationExtent = 0.25*(range(x)[2]-range(x)[1])
# lowX = range(x)[1] - extrapolationExtent
# highX = range(x)[2] + extrapolationExtent
# xPostPred = seq( lowX , highX , length=20 )
# # Define matrix for recording posterior predicted y values at each x value.
# # One row per x value, with each row holding random predicted y values.
# yPostPred = matrix( 0 , nrow=length(xPostPred) , ncol=chainLength )
# # Define matrix for recording HDI limits of posterior predicted y values:
# yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
# # Generate posterior predicted y values.
# # This gets only one y value, at each x, for each step in the chain.
# for ( chainIdx in 1:chainLength ) {
#   yPostPred[,chainIdx] = rnorm( length(xPostPred) ,
#                                 mean = ( mcmcChain[chainIdx,"beta0"] 
#                                          + mcmcChain[chainIdx,"beta1"] * xPostPred 
#                                          + mcmcChain[chainIdx,"beta2"] * xPostPred^2
#                                          # CHANGED
#                                 ) ,
#                                 sd = rep( sigma[chainIdx] , length(xPostPred) ) )
}
