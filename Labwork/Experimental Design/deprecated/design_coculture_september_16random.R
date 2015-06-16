library(vegan)

niveaus <- 10**(c(8, 6, 4))

simulations <- matrix(0, nrow = 100000, ncol = 16)
evenness <- matrix(0, nrow = 100000, ncol = 1)
celaantal <- matrix(0, nrow = 100000, ncol = 1)

#59 # 62 #79 84
set.seed(62)
for(i in 1:100000){
	datapoint <- sample(niveaus, 16, replace = TRUE)
	simulations[i,] <- datapoint
	verhouding <- datapoint/sum(datapoint)
	evenness[i,] <- -sum(verhouding*log(verhouding, base = 16))
	celaantal[i,] <- sum(datapoint)
}



df <- as.data.frame(cbind(simulations, evenness, celaantal))
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12','b13','b14','b15','b16','evenness','celaantal')
# Filter op celaantal
df <- df[log(df$celaantal, base = 10) > 8,] 
# df <- df[log(df$celaantal, base = 10) > 8.1,] # Filter op celaantal


df<-unique(df)
dim(df)

hist(df$evenness)

hist(log(df$celaantal, base = 10))

design <- as.data.frame(matrix(0, nrow = 200, ncol = 18))

#filter op celaantal

design[1:200,] <- df[1:200,]

design <- as.data.frame(design)
names(design) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11','b12','b13','b14','b15','b16','evenness','celaantal')
hist(design$evenness)
hist(design$celaantal)

heatmap(as.matrix(design[,1:10]))

write.csv(design,'design_coculture_september_16random.csv')
cor(design$celaantal, design$evenness)

hist(design$evenness)