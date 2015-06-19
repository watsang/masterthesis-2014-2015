library(texreg)
library(ggplot2)
library(gridExtra)


setwd("D:/univ/2014-2015/thesis/KERMIT/results february")
bacteria <- read.csv("design_results_febr_counts.csv")
bacteria$LMG7866 <- ifelse(bacteria$LMG7866 > 30, 30, bacteria$LMG7866)
bacteria$LMG3203 <- ifelse(bacteria$LMG3203 > 30, 30, bacteria$LMG3203)
bacteria$LMG2954 <- ifelse(bacteria$LMG2954 > 30, 30, bacteria$LMG2954)
bacteria$LMG7878 <- ifelse(bacteria$LMG7878 > 30, 30, bacteria$LMG7878)
head(bacteria)
names(bacteria)


plot(log10(cellcount)~evenness, bacteria, pch = 20)

LMG7866 <- glm(LMG7866~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG2954 <- glm(LMG2954~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG3203 <- glm(LMG3203~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())
LMG7878 <- glm(LMG7878~ evenness+log10(cellcount)+evenness:log10(cellcount), data = bacteria, family = poisson())

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

