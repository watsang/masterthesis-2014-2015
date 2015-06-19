library("MCMCpack")

design <- read.csv("D:/univ/2014-2015/thesis/KERMIT/results february/design_results_febr_counts.csv")
counts <- design$LMG7866
evenness <- design$evenness
cellcount <- log10(design$cellcount)

treatment
posterior <- MCMCpoisson(counts ~ evenness + cellcount, b0 = -30, thin = 5)
plot(posterior)

m1 <- glm.nb(counts ~ evenness + cellcount,-10)
summary(m1)

summary(posterior)

out<-posterior
xyplot(posterior)
xyplot(out) # traceplots
gelman.plot(out) # should be below 1.05 or 1.1
geweke.plot(out) # should all be between -2 and 2
acfplot(out) # should wander around 0

# Results
plot(out)
summary(out)
HPDinterval(out)

### -----------------------------------------
## Using JAGGS
### -----------------------------------------

#-- Example zero-inflated negative binomial regression
zinb <- read.csv("http://www.ats.ucla.edu/stat/data/fish.csv")
zinb
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})

head(zinb)
ggplot(zinb, aes(count, fill = camper)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(camper ~ ., margins=TRUE, scales="free_y")


library(pscl)
m1 <- zeroinfl(count ~ child + camper | persons,
               data = zinb, dist = "negbin", EM = TRUE)
summary(m1)
summary(model <- zeroinfl(counts ~ evenness + cellcount), dist = "negbin", EM = TRUE)

m0 <- update(m1, . ~ 1)

pchisq(2 * (logLik(m1) - logLik(m0)), df = 3, lower.tail=FALSE)

summary(m2 <- glm.nb(count ~ child + camper, data = zinb))

newdata1 <- expand.grid(0:3, factor(0:1), 1:4)
colnames(newdata1) <- c("child", "camper", "persons")
newdata1$phat <- predict(m1, newdata1)

ggplot(newdata1, aes(x = child, y = phat, colour = factor(persons))) +
  geom_point() +
  geom_line() +
  facet_wrap(~camper) +
  labs(x = "Number of Children", y = "Predicted Fish Caught")

####-----------------------------------------
# Example nr. 2-
#########------------------------------------
## data
data("bioChemists", package = "pscl")
## without inflation
## ("art ~ ." is "art ~ fem + mar + kid5 + phd + ment")
fm_pois <- glm(art ~ ., data = bioChemists, family = poisson)
fm_qpois <- glm(art ~ ., data = bioChemists, family = quasipoisson)
fm_nb <- glm.nb(art ~ ., data = bioChemists)
## with simple inflation (no regressors for zero component)
fm_zip <- zeroinfl(art ~ . | 1, data = bioChemists)
fm_zinb <- zeroinfl(art ~ . | 1, data = bioChemists, dist = "negbin")
## inflation with regressors
## ("art ~ . | ." is "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment")
fm_zip2 <- zeroinfl(art ~ . | ., data = bioChemists)
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
head(bioChemists)
bioChemists$


###########################
## Self experiments
#############################
# Source: http://www.inside-r.org/packages/cran/pscl/docs/zeroinfl
# Changed example
setwd("D:/univ/2014-2015/thesis/KERMIT/results february")
bacteria <- read.csv("design_results_febr_counts.csv")

evenness <- bacteria$evenness
cellcount <- bacteria$cellcount
pathogen <- bacteria$LMG7866

bacteria
df <- data.frame(evenness=bacteria$evenness, cellcount=log10(bacteria$cellcount), LMG7866=pathogen)
m1 <- zeroinfl(LMG7866~.|cellcount, data=df, dist = "negbin")
m1 <- hurdle(LMG7866~.|., data=df, dist = "negbin")
summary(m1)

sum(bacteria$LMG7866 == 0)/length(bacteria$LMG7866)

plot(cellcount, pathogen)
x <- order(evenness)
points(cellcount[x],predict(m1)[x], col = "red")
predict1 <- predict(m1, se=T)
predict1$fit

m1 <- zeroinfl(counts ~ evenness + cellcount,
               dist = "negbin", EM = TRUE)
summary(m1)




mean(counts)

n <- length(pathogen)
modelstring <- "
model
{
  # priors
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  
  # likelihood
  for(i in 1:n)
  {
    n50[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1*elev50[i] + beta2*pow(elev50[i],2)
    # this part is here in order to make nice prediction curves:
    prediction[i] ~ dpois(lambda[i])
  } 
}"


