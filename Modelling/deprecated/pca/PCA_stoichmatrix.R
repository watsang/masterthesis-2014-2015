# PCA of example data

library(devtools)
# install_github("ggbiplot", "vqv")

# Load data
# 85963.7 => Helicobacter pylori J99
# 99287.1 => salmonella
# 431946.3 => Escherichia coli SE15

df_helico <-  read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/85963.7.txt")
df_salmon <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/99287.1.txt")
df_ecoli <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/431946.3.txt")

df_salmon$X[1]

ModelSEED.compounds.db <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/ModelSEED-compounds-db.csv")

find_cp <- function(stoichiometric_matrix){
  compound_names <- c()
  for (compound in 1:nrow(stoichiometric_matrix))
  {
    compound_names <- append(compound_names, toString(subset(ModelSEED.compounds.db, DATABASE == toString(stoichiometric_matrix$X[compound]))$ABBREVIATION) )
  }
  return(compound_names)
}

df_helico$X <- find_cp(df_helico)
df_salmon$X <- find_cp(df_salmon)
df_ecoli$X <- find_cp(df_ecoli)


# If you want to show the names on the biplot, run this snippet of code
# column_names <- c()
# for (col in 1:1271){
#   column_names <- append(column_names, toString(df$X[col]))    
# }
# 
# rownames(df) <- column_names
# df

df_helico$X
df_salmon$X
df_ecoli$X

# Analysis Helicobacter
# Perform PCa with compounds as samples
# --------------------------------------------------
# 1) Each sample = compound


df_helico.pca <- prcomp(df_helico[,2:dim(df_helico)], center = T, scale = T)
df_helico.pca

PC1 <- df_helico.pca$rotation[,1]
PC1_eigenvalue <- df_helico.pca$sdev[1]^2

# Check the explained variance
plot(df_helico.pca, type ="l")

# Plot the data on a biplot
biplot(df_helico.pca)

# Show labels on scatterplot
plot(df_helico.pca$x[,1], df_helico.pca$x[,2], xlab ="PC1", ylab = "PC2")
text( df_helico.pca$x[,1], df_helico.pca$x[,2], df_helico$X, pos= 2 )

# Determine the reactions that score positively on PC1
library(prob)
ModelSEED.reactions.db <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/ModelSEED-reactions-db.csv"
PC1names <- names(PC1[PC1 > 0])

# Sort the projected eigenvectors on PC1 from high to low
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ

# Sort the projected eigenvectors on PC1 from low to high
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ


# Analysis salmonnella
# Perform PCa with compounds as samples
# --------------------------------------------------
# 1) Each sample = compound

df_salmon.pca <- prcomp(df_salmon[,2:dim(df_salmon)], center = T, scale = T)
df_salmon.pca

PC1 <- df_salmon.pca$rotation[,1]
PC1_eigenvalue <- df_salmon.pca$sdev[1]^2

# Check the explained variance
plot(df_salmon.pca, type ="l")

# Plot the data on a biplot
biplot(df_salmon.pca)

# Show labels on scatterplot
plot(df_salmon.pca$x[,1], df_salmon.pca$x[,2], xlab ="PC1", ylab = "PC2")
text( df_salmon.pca$x[,1], df_salmon.pca$x[,2], df_salmon$X, pos= 2 )

# Determine the reactions that score positively on PC1
library(prob)
ModelSEED.reactions.db <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/ModelSEED-reactions-db.csv"
                                PC1names <- names(PC1[PC1 > 0])
                                   

# Sort the projected eigenvectors on PC1 from high to low
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ

# Sort the projected eigenvectors on PC1 from low to high
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ

# Analysis ecoli
# Perform PCa with compounds as samples
# --------------------------------------------------
# 1) Each sample = compound

df_ecoli.pca <- prcomp(df_ecoli[,2:dim(df_ecoli)], center = T, scale = T)
df_ecoli.pca

PC1 <- df_ecoli.pca$rotation[,1]
PC1_eigenvalue <- df_ecoli.pca$sdev[1]^2

# Check the explained variance
plot(df_ecoli.pca, type ="l")

# Plot the data on a biplot
biplot(df_ecoli.pca)

# Show labels on scatterplot
plot(df_ecoli.pca$x[,1], df_ecoli.pca$x[,2], xlab ="PC1", ylab = "PC2")
text( df_ecoli.pca$x[,1], df_ecoli.pca$x[,2], df_ecoli$X, pos= 2 )

# Determine the reactions that score positively on PC1
library(prob)
ModelSEED.reactions.db <- read.csv("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/pca/ModelSEED-reactions-db.csv"
                                   PC1names <- names(PC1[PC1 > 0])

# Sort the projected eigenvectors on PC1 from high to low
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ

# Sort the projected eigenvectors on PC1 from low to high
PC1_sorted <- sort(PC1, decreasing = T)
PC1_sorted

# Check in the reactions which ones are most significant
ModelSEED.reactions.db[order(PC1_sorted),][1:50,]$NAME.EQ

                                   