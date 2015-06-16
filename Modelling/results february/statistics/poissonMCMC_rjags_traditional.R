library(rjags)

design <- read.csv("D:/univ/2014-2015/thesis/KERMIT/results february/design_results_febr_counts.csv")
counts <- design$LMG7866
evenness <- design$evenness
cellcount <- log10(design$cellcount)


### -----------------------------------------
## Using JAGGS
### -----------------------------------------

jags.data <- list(N.cells = length(counts), evenness = design$evenness, cellcount = log10(design$cellcount),
                  counts = counts)
jags.data

n <- length(pathogen)
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
    counts[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*evenness[i] + beta2*cellcount[i]
    # this part is here in order to make nice prediction curves:
    prediction[i] ~ dpois(lambda[i])
  } 
}"

cat(modelstring,
    file = "model.txt")

params <- c("beta0", "beta1", "beta2", "prediction")
jm <- jags.model(textConnection(modelstring), 
                 data = jags.data, 
                 n.chains = 3, n.adapt = 1000)

## Run your MCMC iterations!

update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 1)


out <- coda.samples(jm,
                    variable.names=c('counts'),
                    n.iter=1000)
library(lattice)
xyplot(out)
gelman.plot(out) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0
plot(as.mcmc.list(jm.sample$beta0), main = "Beta_0")
plot(as.mcmc.list(jm.sample$beta1), main = "Beta_1")
plot(as.mcmc.list(jm.sample$beta2), main = "Beta_2")
summary(as.mcmc.list(jm.sample$beta2))

#chunk22
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

predict(loess(prds[, 2]~prds[, 1]))

x <- prds[, 1]
y <- prds[, 2]
lo <- loess(prds[, 2]~prds[, 1])
plot(x,y)
lines(x, predict(lo), col='red', lwd=2)


# # # 
# EXTRA CODE
# # #

# treatment
# posterior <- MCMCpoisson(counts ~ evenness + cellcount, b0 = -30, thin = 5)
# plot(posterior)
# 
# m1 <- glm.nb(counts ~ evenness + cellcount,-10)
# summary(m1)
# 
# summary(posterior)
# 
# out<-posterior
# xyplot(posterior)
# xyplot(out) # traceplots
# gelman.plot(out) # should be below 1.05 or 1.1
# geweke.plot(out) # should all be between -2 and 2
# acfplot(out) # should wander around 0
# 
# # Results
# plot(out)
# summary(out)
# HPDinterval(out)