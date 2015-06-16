setwd("~/univ/2014-2015/thesis/initial table3/Re_ preparation co-culture")

library("ineq")

proportions <- read.csv("~/univ/2014-2015/thesis/initial table3/Re_ preparation co-culture/proportions.csv")

gini_prop <- t(apply(proportions, 1, sort))
gini_cum <- t(apply(gini, 1, cumsum))
gini

plot(seq(0,1,.1), seq(0, 1, .1), type = "l")

for (i in 1:194)
{
  lines(Lc(gini_prop[i,]), lwd = 2)
  vector[i]  <- ineq(gini_prop[i,], type="Gini")
}
vector
vector <- rep(NA, 194)
sort(vector)


write.csv(gini, "gini.csv")
