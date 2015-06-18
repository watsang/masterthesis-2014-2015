gini_Wittebolle <- function(vector)
{
  N <- length(vector)
  mu <- mean(vector)
  
  yi <- vector[order(vector, decreasing = T)]
  i <- 1:N
  
  gini <- (N+1)/N - 2 / mu / N^2 * sum(yi *i)

  return(gini)
}

x <- c(1, 1, 1, 1, 1)
x <- rep(100, 1)
x <- sample(niveaus, 10, replace = TRUE)
x
gini_Wittebolle(x)

gini_Thas<-function(x,graph=F) {
  x <- x[order(x, decreasing = T)]
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


niveaus <- c(.4*10**(c(7, 6, 5, 4)), .8*10**(c(7, 6, 5, 4)))

simulations <- matrix(0, nrow = 1000000, ncol = 10)
simulations <- as.data.frame(simulations)
evenness <- matrix(0, nrow = 1000000, ncol = 1)
celaantal <- matrix(0, nrow = 1000000, ncol = 1)

gini_thas <- rep(NA, 10^6)
gini_wittebolle <- rep(NA, 10^6)
set.seed(100)
for(i in 1:1000000){
	datapoint <- sample(niveaus, 10, replace = TRUE)
	simulations[i,] <- datapoint
	verhouding <- datapoint/sum(datapoint)
	evenness[i,] <- -sum(verhouding*log(verhouding, base = 10))
	celaantal[i,] <- sum(datapoint)
  gini_thas[i] <- gini_Thas(simulations[i,])
	gini_wittebolle[i] <- gini_Wittebolle(simulations[i,])	
}


df <- as.data.frame(cbind(simulations, evenness, celaantal, gini_thas, gini_wittebolle))
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal', 'gini_thas', 'gini_wittebolle')
head(df)
par(mfrow = c(1,2), cex = 1.5)
plot(df$evenness, df$celaantal, xlab="Pielou's Evenness", ylab="Cellcount")
rect(0, 10^5.5, 1, 10^7, border = "red", density = NULL, lwd = 5)

plot(df$evenness, df$gini_thas, xlab = "Pielou's evenness", ylab = "Gini_thas")
plot(df$gini_thas, df$celaantal, xlab = "Gini", ylab  = "celaantal")
plot(df$evenness, df$gini_wittebolle, xlab = "Pielou's evenness", ylab = "gini_wittebolle")
plot(df$gini_wittebolle, df$celaantal, xlab = "Gini", ylab  = "celaantal")

df1 <- df[log(df$celaantal, base = 10) < 7,] 
df1 <- df1[log(df$celaantal, base = 10) > 5.5,]
 
df1 <- unique(df1)

# par(mfrow = c(1,1), cex = 1.5)
plot(df1$evenness, df1$celaantal, xlab="evenness", ylab="cellcount")

concentrations <- 10^ceiling(log(df[,1:10], base = 10))

dim(df)

#########################
# Determine volume here #
#########################
# calculate the used volumes

string <- function(x,y=1,z=1) 
{
	as.numeric(substring(toString(x),y,z))*100
}

volumes <- as.data.frame(lapply(df[,1:10],FUN = function(x) {sapply(x,FUN=string)}))

# rowSums(volumes)
vol <- volumes/rowSums(volumes)*1500
#head(df[,1:10]/rowSums(volumes)*1500)
#head(vol)

df[,1:10] <- df[,1:10]/rowSums(volumes)*1500

df <- cbind(df, vol, concentrations)
df$celaantal <- rowSums(df[,1:10])

dim(df)
dim(vol)
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal','gini_thas', 'gini_wittebolle','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
# decorrelate
par(mfrow = c(1,2))
cor(df$evenness, df$celaantal)
plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")
cor(df$gini_thas, df$celaantal)
plot(df$gini_thas, df$celaantal, xlab = "Gini_thas", ylab  = "celaantal")
cor(df$gini_wittebolle, df$celaantal)
plot(df$gini_wittebolle, df$celaantal, xlab = "gini_wittebolle", ylab  = "celaantal")

df <- df[log(df$celaantal, base = 10) < 6.4,] 

hist(df$gini_thas)
hist(df$gini_wittebolle)
hist(df$evenness)

hist(log(df$celaantal, base = 10))

design <- as.data.frame(matrix(0, nrow =194, ncol = 34))

#filter op celaantal

# add 94 random samples
design[101:194,] <- df[1:94,]

hist(df[1:94,]$evenness, main = "", xlab = "", ylab ="")
plot(df[1:94,]$evenness, df[1:94,]$celaantal, main = "", xlab = "", ylab ="")
plot(df[1:94,]$gini_thas, df[1:94,]$celaantal, main = "", xlab = "", ylab ="")
plot(df[1:94,]$gini_wittebolle, df[1:94,]$celaantal, main = "", xlab = "", ylab ="")

# add 100 samples stratified

# 16 samples: high celcount and low evenness
df.sub <- df[log(df$celaantal, base = 10) > 6.2,]
design[1:16,] <- df.sub[order(df.sub$evenness)[1:16],] 

# 8 samples: medium celcount and low evenness
df.sub <- df[(log(df$celaantal, base = 10) < 6.2) & (log(df$celaantal, base = 10) > 5.9),]
design[17:24,] <- df.sub[order(df.sub$evenness)[1:8],] 

# 16 samples: low celcount and low evenness
df.sub <- df[(log(df$celaantal, base = 10) < 5.9),]
design[25:40,] <- df.sub[order(df.sub$evenness)[1:16],] 

# 8 samples: high celcount and medium evenness
df.sub <- df[ (df$evenness < 0.55) & (df$evenness > 0.45),]
df.sub <- df.sub[ log(df.sub$celaantal, base = 10) > 6.2,]
design[41:48,] <- df.sub[1:8,] 

# 4 samples: medium celcount and medium evenness
df.sub <- df[ (df$evenness < 0.55) & (df$evenness > 0.45),]
df.sub <- df.sub[ log(df.sub$celaantal, base = 10) < 6.2,]
df.sub <- df.sub[ log(df.sub$celaantal, base = 10) > 5.9,]
df.sub
design[49:52,] <- df.sub[1:4,]

# 8 samples: low celcount and medium evenness
df.sub <- df[ (df$evenness < 0.55) & (df$evenness > 0.45),]
plot(df.sub$evenness, df.sub$celaantal)
df.sub <- df.sub[ log(df.sub$celaantal, base = 10) < 5.9,]
df.sub
design[53:60,] <- df.sub[1:8,]

# 16 samples: high celcount and high evenness
df.sub <- df[log(df$celaantal, base = 10) > 6.2, ]
design[61:76,] <- df.sub[order(df.sub$evenness, decreasing = T)[1:16],]

# 8 samples: medium celcount and high evenness
df.sub <- df[(log(df$celaantal, base = 10) < 6.2) & (log(df$celaantal, base = 10) > 5.9) , ]
design[77:84,] <- df.sub[order(df.sub$evenness, decreasing = T)[1:8],] 

# 16 samples: low celcount and high evenness
df.sub <- df[(log(df$celaantal, base = 10)  < 5.9),]
design[85:100,] <- df.sub[order(df.sub$evenness, decreasing = T)[1:16],] 


names(design) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal', 'gini_thas', 'gini_wittebolle', 'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
hist(design$evenness, xlab= "evenness", main = "final design")
hist(design[1:100,]$evenness, xlab = "evenness", main = "100 stratified data points")
hist(design[101:194,]$evenness, xlab = "evenness", main = "94 random data points")

hist(design$gini_thas)
hist(design$gini_wittebolle)
hist(design$celaantal)
plot(design$evenness, design$celaantal, xlab="evenness", ylab="cellcount")
plot(design$gini_index, design$celaantal, xlab="gini_index", ylab="cellcount")
plot(design$gini_thas, design$evenness, xlab ="gini_thas", ylab = "pielou's evenness")
plot(design$gini_wittebolle, design$evenness, xlab ="gini_wittebolle", ylab = "pielou's evenness")

par(mfrow = c(1,1))

grid(nx=3,ny=3)
cor(design$gini_index, design$evenness)
heatmap(as.matrix(design[,1:10]))

cor(design$celaantal, design$evenness)
cor(design$evenness, design$gini_wittebolle)
cor(design$evenness, design$gini_thas)


write.csv(design, 'design.csv')

