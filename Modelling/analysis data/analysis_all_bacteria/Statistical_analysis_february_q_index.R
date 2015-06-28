library(texreg)
library(ggplot2)
library(gridExtra)
library(MASS)
install.packages("sqldf")
install.packages("forecast")

setwd("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/q-indices")
bacteria <- read.csv("D:/univ/2014-2015/thesis/KERMIT/results february/design_results_febr_counts.csv")
bacteria$LMG7866 <- ifelse(bacteria$LMG7866 > 30, 30, bacteria$LMG7866)
bacteria$LMG3203 <- ifelse(bacteria$LMG3203 > 30, 30, bacteria$LMG3203)
bacteria$LMG2954 <- ifelse(bacteria$LMG2954 > 30, 30, bacteria$LMG2954)
bacteria$LMG7878 <- ifelse(bacteria$LMG7878 > 30, 30, bacteria$LMG7878)
head(bacteria)
names(bacteria)

q_CCA <- read.table("q_CCA.csv", header=T, sep=",")
q_jacard_compounds <- read.table("q_jacard_compounds.csv", header=T, sep=",")
q_jacard_reactions <- read.table("q_jacard_reactions.csv", header=T, sep=",")
q_jacard_location <- read.table("q_jacard_location.csv", header=T, sep=",")
q_unity <- read.table("q_unity.csv", header=T, sep=",")

# Start calculating residuals for each q_index
calc_residuals <- function(bacteria, similarity_matrix)
{
  residuals <- rep(NA, dim(similarity_matrix)[2])
  for (i in 1:dim(similarity_matrix)[2])
  {
    model1 <- glm(bacteria$LMG7866~ similarity_matrix[,i]+log10(bacteria$cellcount), family = poisson())
    residuals[i] <- sum(abs(model1$residuals))
  }
  return(residuals)
}

model1 <- glm.nb(LMG7866~ evenness+log10(cellcount), data = bacteria)
residuals_nb <- sum(abs(model1$residuals))

Identity_matrix <- calc_residuals(bacteria, q_unity)
CCA <- calc_residuals(bacteria, q_CCA)
Jaccard_compounds <- calc_residuals(bacteria, q_jacard_compounds)
Jaccard_reactions <- calc_residuals(bacteria, q_jacard_reactions)
Jaccard_location <- calc_residuals(bacteria, q_jacard_location)

residuals <- data.frame(q=1:499/100, CCA, Jaccard_compounds, Jaccard_location, 
                        Jaccard_reactions, Identity_matrix, Neg.bin.regression=rep(261, 499))

library(plyr)
library(ggplot2)

library(reshape2)

resid_melt <- melt(residuals, id.vars = 'q')
p_residuals <- ggplot(resid_melt, aes(x=q, y=value, group=variable, col=variable)) + geom_point() +
  theme_grey(base_size = 20, base_family = "") + ylab("Residuals") + 
  guides(colour = guide_legend(override.aes = list(size=4))) + theme(text = element_text(size=25) )
ggsave(file="plot_q_residuals.jpeg")

# ---------------------------------------------------------------------

model2 <- glm(bacteria$LMG7866~ q_unity[,97]+log10(bacteria$cellcount), family = poisson())

plot(bacteria$evenness, predict(model1))
model2 <- glm(bacteria$LMG7866~ q_unity[,97]+log10(bacteria$cellcount), family = poisson())
points(bacteria$evenness, predict(model2), col="red")

plot(evenness, pathogen)
points(bacteria$evenness, predict(model2), col="red")

library(MASS)
m1 <- glm.nb(pathogen ~ evenness + cellcount)
summary(m1)
plot(pathogen~evenness)
points(evenness[order(evenness)], predict(m1)[order(evenness)], col="red")

# rename
names(melted) <- c('id', 'x', 'variable', 'y')

plot(log10(cellcount)~evenness, bacteria, pch = 20)

which.min(residuals)
names(q_unity)[97]

LMG7866 <- glm(LMG7866~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG2954 <- glm(LMG2954~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG3203 <- glm(LMG3203~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG7878 <- glm(LMG7878~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())

# Read csv tables
similarity_CCA <- read.table(file="D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/analysis 10 bacteria/10_bacteria_CCA.csv", header=T, sep = ",", row.names=1)
similarity_CCA <- similarity_CCA / max(max(similarity_CCA))
similarity_CCA <- (similarity_CCA - min(similarity_CCA)) / (max(similarity_CCA) - min(similarity_CCA))
similarity_CCA

head(bacteria)
bacteria_names <- c("pseudo", "bacillus", "serr", "burkho", "para", "entero", "agro", "rhizo", "delftia", "aero")

names(bacteria)[2:11] <- sort(bacteria_names)
bacteria[,2:11] <- bacteria[,2:11][,order(bacteria_names)]
head(bacteria)


# m1 <- extract.glm(model1)
# e
texreg(list(LMG7866, LMG2954, LMG3203, LMG7878), dcolumn = F, booktabs = T, use.packages = FALSE, 
       label = "tab:3", caption = "Two linear models.", float.pos = "h", 
       custom.model.names = c("LMG7866", "LMG2954", "LMG3203", "LMG7878"), scalebox = .5,
       bold = 0.05, caption.above = T)

summary(LMG7866)
summary(LMG2954)
summary(LMG3203)
summary(LMG7878)

model1 <- glm(LMG7866~ evenness+log10(cellcount), data = bacteria, family = poisson())
model2 <- glm(LMG2954~ evenness+log10(cellcount), data = bacteria, family = poisson())
model3 <- glm(LMG3203~ evenness+log10(cellcount), data = bacteria, family = poisson())

texreg(list(model1, model2, model3), dcolumn = F, booktabs = T, use.packages = FALSE, label = "tab:3", caption = "Two linear models.",
       float.pos = "h", custom.model.names = c("LMG7866", "LMG2954", "LMG3203"), 
                                               scalebox = .7, caption.above = T)

summary(model1)
summary(model2)
summary(model3)

sum(abs(model1$residuals))

confint(model1)

par(mfrow = c(2,2))
plot(model1)
plot(model2)

plot(model3)
plot(model4)

plot(bacteria)

aov1 <- aov(LMG7866 ~ replica, data=bacteria) 
aov2 <- aov(LMG3203 ~ replica, data=bacteria) 
aov3 <- aov(LMG2954 ~ replica, data=bacteria) 
aov4 <- aov(LMG7878 ~ replica, data=bacteria) 
summary(aov1)
summary(aov2)
summary(aov3)
summary(aov4)

signif(cor(bacteria)[16:19,], digits = 2)
signif(cor(bacteria)[2:11,16:19], digits = 2)

p1 <- qplot(evenness, LMG7866, data = bacteria, xlab = "Pielou's evenness", main = "LMG7866")
p2 <- qplot(evenness, LMG3203, data = bacteria, xlab = "Pielou's evenness", main = "LMG3203")
p3 <- qplot(evenness, LMG2954, data = bacteria, xlab = "Pielou's evenness", main = "LMG2954")
p4 <- qplot(evenness, LMG7878, data = bacteria, xlab = "Pielou's evenness", main = "LMG7878")
grid.arrange(p1, p2, p3, p4)

fit1 <- predict(model1, type = "response")
fit2 <- predict(model2, type = "response")
plot(bacteria$evenness, fit2)

df1 <- data.frame(evenness = bacteria$evenness, predict = fit1)
df2 <- data.frame(evenness = bacteria$evenness, predict = fit2)

p1 <- ggplot() + 
  geom_point(data=bacteria, aes(evenness, LMG7866), color="black") + 
  geom_point(data=df1, aes(evenness, predict, color="Prediction")) +
  xlab("Pielou's evenness") + ggtitle("LMG7866")

p2 <- ggplot() + 
  geom_point(data=bacteria, aes(evenness, LMG2954), color="black") + 
  geom_point(data=df2, aes(evenness, predict, color="Prediction")) +
  xlab("Pielou's evenness") + ggtitle("LMG2954")
grid.arrange(p1, p2)

