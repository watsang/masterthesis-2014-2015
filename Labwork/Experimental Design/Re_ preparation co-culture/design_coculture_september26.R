gini <- function(x, unbiased = TRUE, na.rm = FALSE){
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind <- is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x <- x[!na.ind]
  n <- length(x)
  mu <- mean(x)
  N <- if (unbiased) n * (n - 1) else n * n
  ox <- x[order(x)]
  dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
  dsum / (mu * N)
}

library(reldist)
x = c(rep(0,9),1)
gini(x, weights = rep(mean(x), length = length(x)))

gini(c(x,1))
gini(x)

niveaus <- c(.4*10**(c(7, 6, 5, 4)), .8*10**(c(7, 6, 5, 4)))

simulations <- matrix(0, nrow = 1000000, ncol = 10)
evenness <- matrix(0, nrow = 1000000, ncol = 1)
celaantal <- matrix(0, nrow = 1000000, ncol = 1)

gini_index <- rep(NA, 10^6)
set.seed(100)
for(i in 1:1000000){
	datapoint <- sample(niveaus, 10, replace = TRUE)
	simulations[i,] <- datapoint
	verhouding <- datapoint/sum(datapoint)
	evenness[i,] <- -sum(verhouding*log(verhouding, base = 10))
	celaantal[i,] <- sum(datapoint)
  gini_index[i] <- ineq(simulations[i,], type = "Gini")
	
}

df <- as.data.frame(cbind(simulations, evenness, celaantal, gini_index))
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal', 'gini_index')
plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")
plot(df$evenness, df$gini_index, xlab = "Pielou's evenness", ylab = "Gini's index")
plot(df$gini_index, df$celaantal, xlab = "Gini", ylab  = "celaantal")

df <- df[log(df$celaantal, base = 10) < 7,] 
df <- df[log(df$celaantal, base = 10) > 5.5,]
 
df <- unique(df)

plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")

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
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal','gini_index','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
# decorrelate
par(mfrow = c(1,2))
cor(df$evenness, df$celaantal)
plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")
cor(df$gini_index, df$celaantal)
plot(df$gini_index, df$celaantal, xlab = "Gini", ylab  = "celaantal")

df <- df[log(df$celaantal, base = 10) < 6.4,] 

hist(gini_index)
hist(df$evenness)

hist(log(df$celaantal, base = 10))

design <- as.data.frame(matrix(0, nrow =194, ncol = 33))

#filter op celaantal

# add 94 random samples
design[101:194,] <- df[1:94,]

hist(df[1:94,]$evenness, main = "", xlab = "", ylab ="")
plot(df[1:94,]$evenness, df[1:94,]$celaantal, main = "", xlab = "", ylab ="")
plot(df[1:94,]$gini_index, df[1:94,]$celaantal, main = "", xlab = "", ylab ="")

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


names(design) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal', 'gini_index', 'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
hist(design$evenness)
hist(design$gini_index)
hist(design$celaantal)
plot(design$evenness, design$celaantal, xlab="evenness", ylab="cellcount")
plot(design$gini_index, design$celaantal, xlab="gini_index", ylab="cellcount")
plot(design$gini_index, design$evenness, xlab ="gini coefficient", ylab = "pielou's evenness")
grid(nx=3,ny=3)
cor(design$gini_index, design$evenness)
heatmap(as.matrix(design[,1:10]))

cor(design$celaantal, design$evenness)
cor(design$celaantal, design$gini_index)

write.csv(design, 'design.csv')

