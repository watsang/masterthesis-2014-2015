niveaus <- c(.4*10**(c(7, 6, 5, 4)), .8*10**(c(7, 6, 5, 4)))

simulations <- matrix(0, nrow = 1000000, ncol = 10)
evenness <- matrix(0, nrow = 1000000, ncol = 1)
celaantal <- matrix(0, nrow = 1000000, ncol = 1)

set.seed(100)
for(i in 1:1000000){
	datapoint <- sample(niveaus, 10, replace = TRUE)
	simulations[i,] <- datapoint
	verhouding <- datapoint/sum(datapoint)
	evenness[i,] <- -sum(verhouding*log(verhouding, base = 10))
	celaantal[i,] <- sum(datapoint)
	
}

df <- as.data.frame(cbind(simulations, evenness, celaantal))
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal')
plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")

hist(c(df[evenness<0.04,][2,1:10]))
typeof(df[evenness<0.04,][2,1:10])

myList <- df[evenness<0.04,][2,1:10]
culture_low_evenness <- c(do.call("cbind",myList))

df1 <- cbind(c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'), culture_low_evenness )
names(df1) <- c("bacteria", "cell count")
df1

hist(df1)
df1 <- df[evenness<0.04,]
write.csv(df1,'df1.csv')

library(reshape2)
new.df<-melt(df1,id.vars="evenness")
names(new.df)=c("evenness","bacteria","rank")

library(ggplot2)
ggplot(new.df, aes(x=evenness,y=rank,fill=bacteria)+
   geom_histogram(stat="identity",position="dodge")

df <- df[log(df$celaantal, base = 10) < 7,] 
df <- df[log(df$celaantal, base = 10) > 5.5,]
 
df <- unique(df)

concentrations <- 10^ceiling(log(df[,1:10], base = 10))

dim(df)

#########################
# Determine volume here #
#########################
# calculate the used volumes

string<-function(x,y=1,z=1) {substring(toString(x),y ,z)}

string <- function(x,y=1,z=1) 
{
	as.numeric(substring(toString(x),y,z))*100
}

test <- as.data.frame(lapply(df[,1:10],FUN = function(x) {sapply(x,FUN=string)}))
test
rowSums(test)
vol <- test/rowSums(test)*1500
head(df[,1:10]/rowSums(test)*1500)

df[,1:10] <- df[,1:10]/rowSums(test)*1500
df <- cbind(df, vol, concentrations)
df$celaantal <- rowSums(df[,1:10])

dim(df)
dim(vol)
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
# decorrelate
cor(df$evenness, df$celaantal)
plot(df$evenness, df$celaantal, xlab="evenness", ylab="celaantal")

df <- df[log(df$celaantal, base = 10) < 6.4,] 

hist(df$evenness)

hist(log(df$celaantal, base = 10))

design <- as.data.frame(matrix(0, nrow =194, ncol = 32))

#filter op celaantal

# add 100 random samples
design[97:194,] <- df[1:98,]

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

design <- as.data.frame(design)
names(design) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10')
hist(design$evenness)
hist(design$celaantal)
plot(design$evenness, design$celaantal, xlab="evenness", ylab="cellcount")
grid(nx=3,ny=3)
heatmap(as.matrix(design[,1:10]))

cor(design$celaantal, design$evenness)

write.csv(design, 'design.csv')

