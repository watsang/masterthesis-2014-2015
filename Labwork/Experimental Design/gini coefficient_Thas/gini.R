gini<-function(x,graph=F) {
  n.sp<-length(x)
  s<-sum(x)
  z0<-1:n.sp*(s/n.sp)
  z1<-cumsum(as.numeric(x))
  if(graph) {
    plot(z0,z1,xlim=c(0,s),ylim=c(0,s))
    abline(c(0,1))
  }
  g<-mean(z1-z0)
  g<-g/mean(z0)
  return(g)
}

vector <- sample(500:1900,1000, replace = T)

x <- c( vector[order(vector, decreasing = T)])
gini(x, graph = T)
A<- rep(3, 8)
B <- cumsum(A)
mean(A-B)
mean(A)