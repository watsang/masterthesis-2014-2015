library(MASS)
library(lmtest)
library(AER)
library(ggplot2)
library(reshape2)
require(foreign)


# setwd("D:/univ/2014-2015/thesis/KERMIT/results february")
setwd("D:/Github/masterthesis-2014-2015/Labwork/Statistical Analysis")
bacteria <- read.csv("design_results_febr_counts.csv")

LMG7866 <- bacteria$LMG7866
LMG7878 <- bacteria$LMG7878
LMG3203 <- bacteria$LMG3203
LMG2954 <- bacteria$LMG2954
evenness <- bacteria$evenness
cellcount <- log10(bacteria$cellcount)
df <- data.frame(evenness, cellcount, LMG7866, LMG7878, LMG3203, LMG2954)

# define loes function
lo <- function(x, y){
  lo <- loess(y~x)
  return(predict(lo, type="response")[order(x)])
}

df2 <- melt(df, id.vars=c("evenness", "cellcount"))
head(df2)
ggplot(df2, aes(x=value, fill=variable)) + geom_histogram(binwidth=1.25) + 
 facet_grid(variable~.) + xlab("Pathogen count") + ylab("Frequency") + 
  theme(axis.text.x = element_text(hjust = 1), axis.text.y = element_text(hjust = 1.5)) +
  theme_set(theme_grey(base_size = 20)) 


test_df <- data.frame(LMG7866=bacteria$LMG7866, LMG3203=bacteria$LMG3203, LMG2954=bacteria$LMG2954, 
                      LMG7878=bacteria$LMG7878, evenness=bacteria$evenness, cellcount=bacteria$cellcount)
test_df <- melt(test_df, id.vars=c("evenness", "cellcount"))
#ggplot(test_df, aes(evenness, value, fill=variable)) + geom_point() + 
 # facet_grid(variable~.)
ggplot(df, aes(x=LMG7866)) + geom_histogram(binwidth=1) + xlab("Pathogen count") + ylab("Frequency") + 
  theme_set(theme_grey(base_size = 30)) 
  
  df

#############################################################################################################
#############################################################################################################
# C L A S S I C     A P P R O A C H
#############################################################################################################
#############################################################################################################
## Poisson
model1 <- glm(LMG7866~evenness+log10(cellcount) +evenness:log10(cellcount),data=df, family = poisson())
model2 <- glm(LMG7878~evenness+log10(cellcount) +evenness:log10(cellcount),data=df, family = poisson())
model3 <- glm(LMG3203~evenness+log10(cellcount) +evenness:log10(cellcount),data=df, family = poisson())
model4 <- glm(LMG2954~evenness+log10(cellcount) +evenness:log10(cellcount),data=df, family = poisson())
summary(model1)
summary(model2)
summary(model3)
summary(model4)

# LMG7878 interaction factor, for the remaining three pathogen the interaction factor is removed

model1 <- glm(LMG7866~ evenness+log10(cellcount), data = bacteria, family = poisson())
model2 <- glm(LMG2954~ evenness+log10(cellcount), data = bacteria, family = poisson())
model3 <- glm(LMG3203~ evenness+log10(cellcount), data = bacteria, family = poisson())

# => C H E C K  F O R  O V E R D I S P E R S I O N
dispersiontest(model1, trafo=2) # overdispersed...
dispersiontest(model2, trafo=2) # overdispersed...
dispersiontest(model3, trafo=2) # overdispersed...

# Source: http://www.ats.ucla.edu/stat/r/dae/nbreg.htm
## negative binomial => used for overdispersed data i.e. variance > mean
model1 <- glm.nb(LMG7866~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria)
model2 <- glm.nb(LMG2954~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria)
model3 <- glm.nb(LMG3203~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria)
summary(model1)
summary(model2)
summary(model3)
# remove interaction factors
model1 <- glm.nb(LMG7866~ evenness+log10(cellcount), data = bacteria)
model2 <- glm.nb(LMG2954~ evenness+log10(cellcount), data = bacteria)
model3 <- glm.nb(LMG3203~ evenness+log10(cellcount), data = bacteria)
summary(model1)
summary(model2)
summary(model3)

texreg(list(model1, model2, model3), dcolumn = F, booktabs = T, use.packages = FALSE, label = "tab:3", caption = "Two linear models.",
       float.pos = "h", custom.model.names = c("LMG7866", "LMG2954", "LMG3203"), 
       scalebox = .7, caption.above = T)


# Check model assumptions of homoscedasticity
par(mfrow = c(2,2))
plot(model1)
plot(model2)
plot(model3)




### check model assumptions
pchisq(2 * (logLik(model_nb) - logLik(model1)), df = 1, lower.tail = FALSE)
# => strongly suggests overdispersion
lrtest(model1, model_nb)

(est <- cbind(Estimate = coef(model_nb), confint(model_nb)))
exp(est)
par(mfrow=c(2,2))
plot(model_nb)

## GOODNESS of FIT
resids1<-residuals(model1, type="pearson") 
sum(resids1^2)
1-pchisq(1915.493,385) 

resids2<-residuals(model_nb, type="pearson") 
sum(resids2^2)
1-pchisq(294.1688,385) 

AIC(model1)
AIC(model_nb)
# For negative binomial regression it has a much smaller AIC == better!
# estimate of the information lost when a given model is used 


## Zero-inflated model: poisson
library(pscl)
model_zi_poisson <- zeroinfl(pathogen~.|., data=df, dist = "poisson")
summary(model_zi_poisson)
plot(pathogen~evenness)
preds <- predict(model1, type = "response", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit
lines(evenness[order(evenness)], lo(evenness, fit), col="red")
lines(evenness[order(evenness)], lo(evenness, upr))
lines(evenness[order(evenness)], lo(evenness, lwr))

## Zero-inflated model: negative binomial => used for overdispersed data i.e. variance > mean
model_zi_negbin <- hurdle(pathogen~.|cellcount, data=df, dist = "negbin")
summary(model_zi_negbin)
plot(pathogen~evenness)
points(evenness, predict(model_zi_negbin), col="red")
prd_zi_negbin <- predict(model_zi_negbin, type = "response", se.fit = TRUE)
prd_zi_negbin
plot(pathogen~evenness)
preds <- predict(model1, type = "response", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit
lines(evenness[order(evenness)], lo(evenness, fit), col="red")
lines(evenness[order(evenness)], lo(evenness, upr))
lines(evenness[order(evenness)], lo(evenness, lwr))

#############################################################################################################
#############################################################################################################
# B A Y E S I A N     A P P R O A C H
#############################################################################################################
#############################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@
# P O I S S O N
#@@@@@@@@@@@@@@@@@@@@@@@@@@
# Simple Metropolis-Hasting
library("MCMCpack")
posterior <- MCMCpoisson(pathogen ~ evenness + cellcount, data=df,
                         beta.start = -10, thin = 5, mcmc=3e5, burnin=1e5)
posterior2 <- MCMCpoisson(pathogen ~ evenness + cellcount, data=df,
                         beta.start = -10, thin = 5, mcmc=3e5, burnin=1e5)
print(BayesFactor(posterior, posterior2))

bind <- mcmc.list(posterior, posterior2)

# MCMCpoisson
plot(posterior)

out<-posterior
xyplot(out) # traceplots
gelman.plot(bind) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0
summary(out)
HPDinterval(out)

# rJAGS GIBBS sampler
library(rjags)

jags.data <- list(N.cells = length(pathogen), evenness = evenness, cellcount = log10(cellcount),
                  counts = pathogen)

# Bayesian analysis with JAGS
modelstring <- "
model
{
  # priors
  beta0 ~ dnorm(0,.01)
  beta_evenness ~ dnorm(0,0.001)
  beta_cellcount ~ dnorm(0,.01)
  
  # likelihood
  for(i in 1:N.cells)
  {
    counts[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta_evenness*evenness[i] + beta_cellcount*cellcount[i]
    # this part is here in order to make nice prediction curves:
    prediction[i] ~ dpois(lambda[i])
  } 
}"

# Initiate the model
initial <- list(list(beta0=3,beta_evenness=-10, beta_cellcount=-1),
                list(beta0=10,beta_evenness=-10, beta_cellcount=-10),
                list(beta0=10,beta_evenness=-10, beta_cellcount=-10))

params <- c("beta0", "beta_evenness", "beta_cellcount", "prediction")
jm <- jags.model(textConnection(modelstring), inits=initial,
                 data = jags.data,
                 n.chains = 3, n.adapt = 2000)

# Update initial model 
update(jm, 1000)

# Run MCMC chain
out <- coda.samples(jm,
                    variable.names=params[1:3],
                    n.iter=10000, thin= 1)
out_temp <- out
geweke.diag(out)
# Diagnostics checks
library(lattice)
xyplot(out) # traceplots => BAD MIXING characteristics 
gelman.plot(out) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0
autocorr.plot(out)
rejectionRate(out)

# Results
plot(out)
summary(out)
HPDinterval(out)

#@@@@@@@@@@@@@@@@@@@@@@@@@@
# N E G A T I V E  B I N O M I A L  R E G R E S S I O N
#@@@@@@@@@@@@@@@@@@@@@@@@@@
jags.data <- list(N.cells = length(pathogen), evenness = evenness, cellcount = log10(cellcount),
                 counts = pathogen)

# Bayesian analysis with JAGS
modelstring <- "
model
{
# likelihood
for(i in 1:N.cells){
  mu[i] <- beta0 + beta_evenness*evenness[i] + beta_cellcount*cellcount[i]
  lambda[i] <- exp(mu[i])
  p[i] <- abs(r/(r+lambda[i]))
  # this part is here in order to make nice prediction curves:
  counts[i] ~ dnegbin(p[i],r)
}

# priors
beta0 ~ dnorm(0,.001)
beta_evenness ~ dnorm(0,0.001)
beta_cellcount ~ dnorm(0,.001)
r ~ dunif(0,20)

}"

# Initiate the model
initial <- list(list(beta0=0,beta_evenness=-2, beta_cellcount=0, r =1),
                list(beta0=0,beta_evenness=0, beta_cellcount=0, r =1),
                list(beta0=0,beta_evenness=2, beta_cellcount=0, r =1))
params <- c("beta0", "beta_evenness", "beta_cellcount", "r")
jm <- jags.model(textConnection(modelstring), inits=initial,
                 data = jags.data,
                 n.chains = 3, n.adapt = 10000)

# Update initial model 
update(jm, 5000)

# Run MCMC chain
out <- coda.samples(jm,
                    variable.names=params[1:3],
                    n.iter=50000, thin= 1)
out_temp <- out
geweke.diag(out)
# Diagnostics checks
library(lattice)
xyplot(out) # traceplots => BAD MIXING characteristics 
gelman.plot(out) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0
autocorr.plot(out)
rejectionRate(out)

# Results
plot(out)
summary(out)
HPDinterval(out)

jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 1)

predictions <- summary(as.mcmc.list(jm.sample$prediction))
prds <- data.frame(sc50 = evenness, predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

scale2 <- function(x) {
  sdx <- sqrt(var(x))
  meanx <- mean(x)
  return((x - meanx)/sdx)
}

plot(evenness, counts, cex = 1, col = "lightgrey", pch = 19, ylab = "Colony count", 
     xlab = "Evennes")
x <- prds[, 1]
lo1 <- loess(prds[, 2]~prds[, 1])
lo2 <- loess(prds[, 4]~prds[, 1])
lo3 <- loess(prds[, 6]~prds[, 1])
lines(x, predict(lo1), lwd = 2)
lines(x, predict(lo2),lwd = 2, col = "red")
lines(x, predict(lo3), lwd = 2)


########################################
# T R Y  A G A I N
######################################3


# B i b l i o g r a p h y 
#### ASSUMPTIONS
# 
# Source: http://www.scalelive.com/negative-binomial-regression.html
# Residuals can be thought of as the error associated with predicting or estimating outcomes using predictor
# variables.  Residual analysis is extremely important for meeting the linearity, normality, and homogeneity 
# of variance assumptions of negative binomial regression. 