niveaus <- 10**(c(8, 6, 4))

simulations <- matrix(0, nrow = 100000, ncol = 10)
evenness <- matrix(0, nrow = 100000, ncol = 1)
celaantal <- matrix(0, nrow = 100000, ncol = 1)

set.seed(1)
for(i in 1:100000){
	datapoint <- sample(niveaus, 10, replace = TRUE)
	simulations[i,] <- datapoint
	verhouding <- datapoint/sum(datapoint)
	evenness[i,] <- -sum(verhouding*log(verhouding, base = 10))
	celaantal[i,] <- sum(datapoint)
	
}

df <- as.data.frame(cbind(simulations, evenness, celaantal))
names(df) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal')
# Filter op celaantal
df <- df[log(df$celaantal, base = 10) > 8,] 
# df <- df[log(df$celaantal, base = 10) > 8.1,] # Filter op celaantal


df<-unique(df)
dim(df)

hist(df$evenness)

hist(log(df$celaantal, base = 10))


design <- as.data.frame(matrix(0, nrow = 200, ncol = 12))

#filter op celaantal

design[1:100,] <- df[1:100,]

#build design with 100 samples high and low level + 100 samples randomized

design[101:150,] <- df [order(df$evenness)[1:50],]

design[151:200,] <- df [order(df$evenness),][(length(df$evenness)-49):length(df$evenness),]

design <- as.data.frame(design)
names(design) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal')
hist(design$evenness)
hist(design$celaantal)

heatmap(as.matrix(design[,1:10]))

#build design.stratified1 with 80 samples high and low level + 120 samples stratified around 0.25 0.5 and 0.75

design.stratified1 <- as.data.frame(matrix(0, nrow = 200, ncol = 12))

design.stratified1[1:40,] <- df [order(df$evenness)[1:40],]

design.stratified1[41:80,] <- df [order(df$evenness),][(length(df$evenness)-39):length(df$evenness),]

design.stratified1[81:120,] <- df[df$evenness > 0.2 & df$evenness < 0.34,][1:40,]  # Zoek op evenness = .3

design.stratified1[121:160,] <- df[df$evenness > 0.475 & df$evenness < 0.525,][1:40, ] # Filter op evenness = .5

design.stratified1[161:200,] <- df[df$evenness > 0.7 & df$evenness < 0.8,][1:40, ] # Filter op evenness = .75

design.stratified1 <- as.data.frame(design.stratified1)
names(design.stratified1) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal')
hist(design.stratified1$evenness, breaks=seq(0,1,by=0.05))
hist(design.stratified1$celaantal)

#build design.stratified2 with 80 samples high and low level + 120 samples stratified around 0.2, 0.35, 0.5, 0.65, 0.8

design.stratified2 <- as.data.frame(matrix(0, nrow = 200, ncol = 12))

design.stratified2[1:40,] <- df [order(df$evenness)[1:40],]

design.stratified2[41:80,] <- df [order(df$evenness),][(length(df$evenness)-39):length(df$evenness),]

#build design.stratified2 with 80 samples high and low level + 120 samples stratified around 0.25 0.5 and 0.75

design.stratified2[81:104,] <- df[df$evenness > 0.18 & df$evenness < 0.25,][1:24,]  # Zoek op evenness = .2

design.stratified2[105:128,] <- df[df$evenness > 0.25 & df$evenness < 0.34,][1:24,]  # Zoek op evenness = .35

design.stratified2[129:152,] <- df[df$evenness > 0.485 & df$evenness < 0.515,][1:24, ] # Filter op evenness = .5

design.stratified2[153:176,] <- df[df$evenness > 0.625 & df$evenness < 0.675,][1:24, ] # Filter op evenness = .65

design.stratified2[177:200,] <- df[df$evenness > 0.75 & df$evenness < 0.84,][1:24, ] # Filter op evenness = .8

design.stratified2 <- as.data.frame(design.stratified2)
names(design.stratified2) <- c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','evenness','celaantal')
hist(design.stratified2$evenness, breaks=seq(0,1,by=0.05))
hist(design.stratified2$celaantal)

cor(design$celaantal, design$evenness)
