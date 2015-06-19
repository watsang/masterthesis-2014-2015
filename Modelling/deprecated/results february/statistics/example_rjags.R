library(spatstat)
library(raster)
library(sp)
library(googleVis)
library(rgdal)
library(knitr)
opts_chunk$set(cache = FALSE)

#chunk2
plot(bei$x, bei$y, pch = 19, cex = 0.5, main = "Spatial distribution of individuals in the 50-ha Barro Colorado plot", 
     xlab = "x coordinate [m]", ylab = "y coordinate [m]", frame = FALSE)
abline(h = 0, col = "grey")
abline(h = 500, col = "grey")
abline(v = 0, col = "grey")
abline(v = 1000, col = "grey")

#chunk3
# coarsening the predictor data into the 50 x 50 m grid by taking the mean
# of the 5 x 5 m grid cells:
elev <- raster(bei.extra[[1]])

# cropping the data so that they have exactly 500 x 1000 cells
ext <- extent(2.5, 1002.5, 2.5, 1002.5)
elev <- crop(elev, ext)

# aggregating the elevation data
elev50.raster <- aggregate(elev, fact = 10, fun = mean)

# fitting the point data into the 50 x 50 m grid
xy <- data.frame(x = bei$x, y = bei$y)
n50.raster <- rasterize(xy, elev50.raster, fun = "count")

# replacing the NA values by 0
n50.raster[is.na(n50.raster)] <- 0

#chunk4
plot(stack(elev50.raster, n50.raster), main = c("Predictor: Mean Elevation in 50x50 m cells", 
                                                "Response: # of Individuals in 50x50 m cells"), axes = FALSE)

#chunk5
plot(elev50.raster[], n50.raster[], cex = 1, pch = 19, col = "grey", ylab = "# of Individuals", 
     xlab = "Mean Elevation [m]")


#chunk6
scale2 <- function(x) {
  sdx <- sqrt(var(x))
  meanx <- mean(x)
  return((x - meanx)/sdx)
}

elev50 <- scale2(elev50.raster[])

pow.elev50 <- elev50^2  # (I will be fitting a polynomial)
n50 <- n50.raster[]




#chunk7
explore.model <- function(beta0, beta1, beta2, elev50.raster) {
  elev50 <- scale2(elev50.raster[])
  n = length(elev50)
  ## parameter lambda of the poisson distribution
  lambda <- exp(beta0 + beta1 * elev50 + beta2 * (elev50^2))
  ## 95% prediction intervals calculated by the qpois() function
  lambda.up95 <- qpois(rep(0.975, times = n), lambda)
  lambda.low95 <- qpois(rep(0.025, times = n), lambda)
  ## new data simulated by rpois() function
  samples <- rpois(n = n, lambda = lambda)
  
  ## plotting it all
  par(mfrow = c(1, 2))
  plot(elev50, samples, pch = 19, col = "grey", xlab = "Scaled Mean Elevation", 
       ylab = "# of Individuals")
  points(elev50, lambda, pch = 19, col = "red")
  points(elev50, lambda.up95, pch = 19, col = "red", cex = 0.5)
  points(elev50, lambda.low95, pch = 19, col = "red", cex = 0.5)
  
  elev50.raster[] <- lambda
  plot(elev50.raster, main = "Predicted Lambda (mean # of individuals)", axes = FALSE)
}

#chunk8
explore.model(beta0 = 3, beta1 = 0.5, beta2 = -0.5, elev50.raster = elev50.raster)

#chunk9
m.glm <- glm(n50 ~ elev50 + pow.elev50, family = "poisson")
summary(m.glm)

#chunk10
elev.seq <- seq(-3, 2, by = 0.05)
new.data <- data.frame(elev50 = elev.seq, pow.elev50 = elev.seq^2)
new.predict <- predict(m.glm, newdata = new.data, type = "response")

#chunk11
plot(elev50, n50, cex = 1, col = "lightgrey", pch = 19, ylab = "# of Individuals", 
     xlab = "Scaled Mean Elevation")
lines(elev.seq, new.predict, col = "red", lwd = 2)


#chunk12
LogLike <- function(dat, par) {
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  # the deterministic part of the model:
  lambda <- exp(beta0 + beta1 * dat$x + beta2 * (dat$x^2))
  # and here comes the negative log-likelihood of the whole dataset, given the
  # model:
  LL <- -sum(dpois(dat$y, lambda, log = TRUE))
  return(LL)
}


#chunk13
dat <- data.frame(y = n50, x = elev50)

#chunk14
par <- c(beta0 = 10, beta1 = 0.05, beta2 = -0.5)
LogLike(dat = dat, par = par)


#chunk15
beta0 <- rnorm(1)
beta1 <- rnorm(1)
beta2 <- rnorm(1)
par <- c(beta0, beta1, beta2)

#chunk16
m.like <- optim(par = par, fn = LogLike, dat = dat)
m.like

#chunk11
plot(elev50, n50, cex = 1, col = "lightgrey", pch = 19, ylab = "# of Individuals", 
     xlab = "Scaled Mean Elevation")
new.predict <- exp(m.like$par[1] + m.like$par[2] * elev.seq + m.like$par[3] * 
                     (elev.seq^2))
lines(elev.seq, new.predict, col = "red", lwd = 2)


#chunk17
library(rjags)

#chunk18
jags.data <- list(N.cells = length(n50), n50 = n50, elev50 = elev50)

#chunk19
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
      }
"

#chunk20
params <- c("beta0", "beta1", "beta2", "prediction")
jm <- jags.model("model.txt", data = jags.data, n.chains = 3, n.adapt = 1000)

#chunk21
update(jm, n.iter = 1000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 1000, thin = 1)
plot(as.mcmc.list(jm.sample$beta0), main = "Beta_0")
plot(as.mcmc.list(jm.sample$beta1), main = "Beta_1")
plot(as.mcmc.list(jm.sample$beta2), main = "Beta_2")
summary(as.mcmc.list(jm.sample$beta2))


#chunk22
predictions <- summary(as.mcmc.list(jm.sample$prediction))
prds <- data.frame(sc50 = scale2(elev50), predictions$quantiles)
prds <- prds[order(prds[, 1]), ]

#chunk23
plot(scale2(elev50), n50, cex = 1, col = "lightgrey", pch = 19, ylab = "# of Individuals", 
     xlab = "Scaled Mean Elevation")
lines(prds[, 1], prds[, 2], lwd = 2)
lines(prds[, 1], prds[, 4], lwd = 2, col = "red")
lines(prds[, 1], prds[, 6], lwd = 2)

legend("topleft", legend = c("95% P.I.", "lambda_i"), col = c("black", "red"), 
       lwd = c(2, 2))

