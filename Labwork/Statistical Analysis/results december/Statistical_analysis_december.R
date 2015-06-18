bacteria <- read.csv('Gini-table_converted.csv')
bacteria <- read.csv("~/univ/2014-2015/thesis/results december/design_results.csv")
# bacteria$count <- ifelse(bacteria$count == +, 30, bacteria$count)
typeof(bacteria$count)
bacteria$count
summary(bacteria)

plot(log(celaantal)~evenness, bacteria, pch = 20)

model <- glm(count~ evenness+log(celaantal)+evenness:log(celaantal), data = bacteria, family = poisson())
summary(model)
confint(model)

par(mfrow = c(2,2))
plot(model)