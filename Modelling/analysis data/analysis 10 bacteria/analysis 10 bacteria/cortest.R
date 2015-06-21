setwd("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/analysis 10 bacteria")

data(iris)

iris

cov(iris[,1:4])

b <- iris[,1:4]
for (i in 1:4)
{
	for (j in 1:dim(iris)[1])
	{
		a <- iris[,i] - mean(iris[,i])
		b[j,i] <-  / (sd(iris[,i])))
	}
}

as.matrix(t(b))  %*% as.matrix(b)

b <- iris[,1:4]
b <- b - lapply(iris[,1:4], mean)
b <- b / 

min(b,

lapply(iris[,1:4], max) - lapply(iris[,1:4], min)

b



