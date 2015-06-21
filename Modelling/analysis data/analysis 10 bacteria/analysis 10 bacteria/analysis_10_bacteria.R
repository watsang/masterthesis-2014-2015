

# library("extrafont")
# font_import()
# loadfonts()
# loadfonts(device = "postscript")
# OBJECTIVE
# ==============
#
# 1)  Calculate S^T and S for organisms
# 
# 2)  Perform PCA on both S and S^T
#
# 3) Compare PCAs
#
#

# ===========================================================

working_directory <- "D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining"
setwd(working_directory)
dir.create(file.path(working_directory, "NAME_EQUATION_PCA"), showWarnings = FALSE)
#setwd("~/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/analysis 10 bacteria")
df_rhizobium <- read.csv("data_organism_reactions/Stoichiometric Matrices/394.7.txt")
df_agro <- read.csv("data_organism_reactions/Stoichiometric Matrices/176299.3.txt")
df_burkho <- read.csv("data_organism_reactions/Stoichiometric Matrices/269483.3.txt")
df_para <- read.csv("data_organism_reactions/Stoichiometric Matrices/318586.5.txt")
df_bacillus <- read.csv("data_organism_reactions/Stoichiometric Matrices/388400.4.txt")

# Change Lipid A disaccharide to Lipid A2
ModelSEED.compounds.db <- read.csv("data_preprocessing/ModelSEED-compounds-db.csv")
ModelSEED.reactions.db <- read.csv("data_preprocessing/ModelSEED-reactions-db.csv")

find_cp <- function(stoichiometric_matrix){
  compound_names <- c()
  for (compound in 1:nrow(stoichiometric_matrix))
  {
    compound_names <- c(compound_names, toString(subset(ModelSEED.compounds.db, DATABASE == toString(stoichiometric_matrix$X[compound]) )$ABBREVIATION) )
  }
  return(compound_names)
}

find_reaction <- function(stoichiometric_matrix){
  reaction_names <- c()
  for (reaction in 1:dim(stoichiometric_matrix)[2])
  {
    reaction_names <- c(reaction_names, toString(subset(ModelSEED.reactions.db, DATABASE == toString(names(stoichiometric_matrix)[reaction]))$NAME) )
  }
  return(reaction_names)
}

find_full_reaction <- function(stoichiometric_matrix){
  reaction_names <- c()
  for (reaction in 1:dim(stoichiometric_matrix)[2])
  {
    reaction_names <- c(reaction_names, toString(subset(ModelSEED.reactions.db, DATABASE == toString(names(stoichiometric_matrix)[reaction]))$NAME.EQ) )
  }
  return(reaction_names)
}

df_rhizo.NAME.EQ <- find_full_reaction(df_rhizobium[,2:dim(df_rhizobium)[2]])
df_agro.NAME.EQ <- find_full_reaction(df_agro[,2:dim(df_agro)[2]])
df_para.NAME.EQ <- find_full_reaction(df_para[,2:dim(df_para)[2]])
df_burkho.NAME.EQ <- find_full_reaction(df_burkho[,2:dim(df_burkho)[2]])
df_bacillus.NAME.EQ <- find_full_reaction(df_bacillus[,2:dim(df_bacillus)[2]])

# Switch compound code with real compound name
# Switch reaction code with reaction name


names(df_rhizobium) <- c("X", find_reaction(df_rhizobium[,2:dim(df_rhizobium)[2]]))
names(df_agro) <- c("X", find_reaction(df_agro[,2:dim(df_agro)[2]]))
names(df_para) <- c("X", find_reaction(df_para[,2:dim(df_para)[2]]))
names(df_burkho) <- c("X", find_reaction(df_burkho[,2:dim(df_burkho)[2]]))
names(df_bacillus) <- c("X", find_reaction(df_bacillus[,2:dim(df_bacillus)[2]]))

df_rhizobium$X <- find_cp(df_rhizobium)
df_agro$X <- find_cp(df_agro)
df_burkho$X <- find_cp(df_burkho)
df_para$X <- find_cp(df_para)
df_bacillus$X <- find_cp(df_bacillus)

# Calculate Stoichiometric matrix * stoichiometric matrix ^ T
# and Stoichiometric matrix ^ T * stoichiometric matrix 
# ------------------------------------------------------------
S_rhizo <- data.matrix( df_rhizobium[,2:dim(df_rhizobium)[2]] )
row.names(S_rhizo) <- df_rhizobium$X

S_agro <- data.matrix( df_agro[,2:dim(df_agro)[2]] )
row.names(S_agro) <- df_agro$X

S_burkho <- data.matrix( df_burkho[,2:dim(df_burkho)[2]] )
row.names(S_burkho) <- df_burkho$X

S_para <- data.matrix( df_para[,2:dim(df_para)[2]] )
row.names(S_para) <- df_para$X

S_bacillus <- data.matrix( df_bacillus[,2:dim(df_bacillus)[2]] )
row.names(S_bacillus) <- df_bacillus$X

ST_bacillus <- t(S_bacillus)
ST_rhizo <- t(S_rhizo)
ST_agro <- t(S_agro)
ST_para <- t(S_para)
ST_burkho <- t(S_burkho)

# PCA
# ===============
# First calculate PCA for all organisms, then compare
#
# 1: Calculations
#---------------
# A: Perform PCA on all S
S_burkho.pca <- prcomp(S_burkho, center = T, scale = T)
selection <- colSums(S_burkho) == 0
write.csv(S_burkho, "S_burkho.csv" )
S_rhizo.pca <- prcomp(S_rhizo, center = T, scale = T)
S_agro.pca <- prcomp(S_agro, center = T, scale = T)
S_para.pca <- prcomp(S_para, center = T, scale = T)
S_bacillus.pca <- prcomp(S_bacillus, center = T, scale = T)

# B: select eigenvectors of respective principle components
PC1.S_burkho <- S_burkho.pca$rotation[,1]
PC2.S_burkho <- S_burkho.pca$rotation[,2]
PC3.S_burkho <- S_burkho.pca$rotation[,3]
PC4.S_burkho <- S_burkho.pca$rotation[,4]
PC5.S_burkho <- S_burkho.pca$rotation[,5]
PC6.S_burkho <- S_burkho.pca$rotation[,6]

PC1.eigenvalue.S_burkho <- S_burkho.pca$sdev[1]^2
PC2.eigenvalue.S_burkho <- S_burkho.pca$sdev[2]^2
PC3.eigenvalue.S_burkho <- S_burkho.pca$sdev[3]^2
PC4.eigenvalue.S_burkho <- S_burkho.pca$sdev[4]^2
PC5.eigenvalue.S_burkho <- S_burkho.pca$sdev[5]^2
PC6.eigenvalue.S_burkho <- S_burkho.pca$sdev[6]^2

PC1.S_rhizo <- S_rhizo.pca$rotation[,1]
PC2.S_rhizo <- S_rhizo.pca$rotation[,2]
PC3.S_rhizo <- S_rhizo.pca$rotation[,3]
PC4.S_rhizo <- S_rhizo.pca$rotation[,4]
PC5.S_rhizo <- S_rhizo.pca$rotation[,5]
PC6.S_rhizo <- S_rhizo.pca$rotation[,6]

PC1.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[1]^2
PC2.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[2]^2
PC3.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[3]^2
PC4.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[4]^2
PC5.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[5]^2
PC6.eigenvalue.S_rhizo <- S_rhizo.pca$sdev[6]^2

PC1.S_agro <- S_agro.pca$rotation[,1]
PC2.S_agro <- S_agro.pca$rotation[,2]
PC3.S_agro <- S_agro.pca$rotation[,3]
PC4.S_agro <- S_agro.pca$rotation[,4]
PC5.S_agro <- S_agro.pca$rotation[,5]
PC6.S_agro <- S_agro.pca$rotation[,6]

PC1.eigenvalue.S_agro <- S_agro.pca$sdev[1]^2
PC2.eigenvalue.S_agro <- S_agro.pca$sdev[2]^2
PC3.eigenvalue.S_agro <- S_agro.pca$sdev[3]^2
PC4.eigenvalue.S_agro <- S_agro.pca$sdev[4]^2
PC5.eigenvalue.S_agro <- S_agro.pca$sdev[5]^2
PC6.eigenvalue.S_agro <- S_agro.pca$sdev[6]^2

PC1.S_para <- S_para.pca$rotation[,1]
PC2.S_para <- S_para.pca$rotation[,2]
PC3.S_para <- S_para.pca$rotation[,3]
PC4.S_para <- S_para.pca$rotation[,4]
PC5.S_para <- S_para.pca$rotation[,5]
PC6.S_para <- S_para.pca$rotation[,6]

PC1.eigenvalue.S_para <- S_para.pca$sdev[1]^2
PC2.eigenvalue.S_para <- S_para.pca$sdev[2]^2
PC3.eigenvalue.S_para <- S_para.pca$sdev[3]^2
PC4.eigenvalue.S_para <- S_para.pca$sdev[4]^2
PC5.eigenvalue.S_para <- S_para.pca$sdev[5]^2
PC6.eigenvalue.S_para <- S_para.pca$sdev[6]^2

PC1.S_bacillus <- S_bacillus.pca$rotation[,1]
PC2.S_bacillus <- S_bacillus.pca$rotation[,2]
PC3.S_bacillus <- S_bacillus.pca$rotation[,3]
PC4.S_bacillus <- S_bacillus.pca$rotation[,4]
PC5.S_bacillus <- S_bacillus.pca$rotation[,5]
PC6.S_bacillus <- S_bacillus.pca$rotation[,6]

PC1.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[1]^2
PC2.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[2]^2
PC3.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[3]^2
PC4.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[4]^2
PC5.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[5]^2
PC6.eigenvalue.S_bacillus <- S_bacillus.pca$sdev[6]^2

ST_burkho.pca <- prcomp(ST_burkho, center = T, scale = T)
ST_rhizo.pca <- prcomp(ST_rhizo, center = T, scale = T)
ST_agro.pca <- prcomp(ST_agro, center = T, scale = T)
ST_para.pca <- prcomp(ST_para, center = T, scale = T)
ST_bacillus.pca <- prcomp(ST_bacillus, center = T, scale = T)

# B: select eigenvectors of respective principle components
PC1.ST_burkho <- ST_burkho.pca$rotation[,1]
PC2.ST_burkho <- ST_burkho.pca$rotation[,2]
PC3.ST_burkho <- ST_burkho.pca$rotation[,3]
PC4.ST_burkho <- ST_burkho.pca$rotation[,4]
PC5.ST_burkho <- ST_burkho.pca$rotation[,5]
PC6.ST_burkho <- ST_burkho.pca$rotation[,6]

PC1.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[1]^2
PC2.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[2]^2
PC3.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[3]^2
PC4.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[4]^2
PC5.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[5]^2
PC6.eigenvalue.ST_burkho <- ST_burkho.pca$sdev[6]^2

PC1.ST_rhizo <- ST_rhizo.pca$rotation[,1]
PC2.ST_rhizo <- ST_rhizo.pca$rotation[,2]
PC3.ST_rhizo <- ST_rhizo.pca$rotation[,3]
PC4.ST_rhizo <- ST_rhizo.pca$rotation[,4]
PC5.ST_rhizo <- ST_rhizo.pca$rotation[,5]
PC6.ST_rhizo <- ST_rhizo.pca$rotation[,6]

PC1.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[1]^2
PC2.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[2]^2
PC3.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[3]^2
PC4.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[4]^2
PC5.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[5]^2
PC6.eigenvalue.ST_rhizo <- ST_rhizo.pca$sdev[6]^2

PC1.ST_agro <- ST_agro.pca$rotation[,1]
PC2.ST_agro <- ST_agro.pca$rotation[,2]
PC3.ST_agro <- ST_agro.pca$rotation[,3]
PC4.ST_agro <- ST_agro.pca$rotation[,4]
PC5.ST_agro <- ST_agro.pca$rotation[,5]
PC6.ST_agro <- ST_agro.pca$rotation[,6]

PC1.eigenvalue.ST_agro <- ST_agro.pca$sdev[1]^2
PC2.eigenvalue.ST_agro <- ST_agro.pca$sdev[2]^2
PC3.eigenvalue.ST_agro <- ST_agro.pca$sdev[3]^2
PC4.eigenvalue.ST_agro <- ST_agro.pca$sdev[4]^2
PC5.eigenvalue.ST_agro <- ST_agro.pca$sdev[5]^2
PC6.eigenvalue.ST_agro <- ST_agro.pca$sdev[6]^2

PC1.ST_para <- ST_para.pca$rotation[,1]
PC2.ST_para <- ST_para.pca$rotation[,2]
PC3.ST_para <- ST_para.pca$rotation[,3]
PC4.ST_para <- ST_para.pca$rotation[,4]
PC5.ST_para <- ST_para.pca$rotation[,5]
PC6.ST_para <- ST_para.pca$rotation[,6]

PC1.eigenvalue.ST_para <- ST_para.pca$sdev[1]^2
PC2.eigenvalue.ST_para <- ST_para.pca$sdev[2]^2
PC3.eigenvalue.ST_para <- ST_para.pca$sdev[3]^2
PC4.eigenvalue.ST_para <- ST_para.pca$sdev[4]^2
PC5.eigenvalue.ST_para <- ST_para.pca$sdev[5]^2
PC6.eigenvalue.ST_para <- ST_para.pca$sdev[6]^2

PC1.ST_bacillus <- ST_bacillus.pca$rotation[,1]
PC2.ST_bacillus <- ST_bacillus.pca$rotation[,2]
PC3.ST_bacillus <- ST_bacillus.pca$rotation[,3]
PC4.ST_bacillus <- ST_bacillus.pca$rotation[,4]
PC5.ST_bacillus <- ST_bacillus.pca$rotation[,5]
PC6.ST_bacillus <- ST_bacillus.pca$rotation[,6]

PC1.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[1]^2
PC2.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[2]^2
PC3.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[3]^2
PC4.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[4]^2
PC5.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[5]^2
PC6.eigenvalue.ST_bacillus <- ST_bacillus.pca$sdev[6]^2

# 2: Comparison
# -------------
# Plot PC 2 and 3
library(ggplot2)
library(gridExtra)

## Text selection

farthest_distance <- function(vector1, vector2)
{
  distance <- (vector1^2 + vector2^2)^.5  
  sorted_distance <- order(distance, decreasing = T)
  selection <- sorted_distance[1:15] # Select the first
  return(selection)
}

# Create datasets for rotated data points
rotated_S_burkho <- as.data.frame(S_burkho.pca$x)
rotated_S_agro <- as.data.frame(S_agro.pca$x)
rotated_S_para <- as.data.frame(S_para.pca$x)
rotated_S_bacillus <- as.data.frame(S_bacillus.pca$x)
rotated_S_rhizo <- as.data.frame(S_rhizo.pca$x)


# Make a selection of resp PC's
selection_S_burkho_12 <- farthest_distance(rotated_S_burkho$PC1, rotated_S_burkho$PC2)
selection_S_para_12 <- farthest_distance(rotated_S_para$PC1, rotated_S_para$PC2)
selection_S_agro_12 <- farthest_distance(rotated_S_agro$PC1, rotated_S_agro$PC2)
selection_S_bacillus_12 <- farthest_distance(rotated_S_bacillus$PC1, rotated_S_bacillus$PC2)
selection_S_rhizo_12 <- farthest_distance(rotated_S_rhizo$PC1, rotated_S_rhizo$PC2)

selection_S_burkho_23 <- farthest_distance(rotated_S_burkho$PC2, rotated_S_burkho$PC3)
selection_S_para_23 <- farthest_distance(rotated_S_para$PC2, rotated_S_para$PC3)
selection_S_agro_23 <- farthest_distance(rotated_S_agro$PC2, rotated_S_agro$PC3)
selection_S_bacillus_23 <- farthest_distance(rotated_S_bacillus$PC2, rotated_S_bacillus$PC3)
selection_S_rhizo_23 <- farthest_distance(rotated_S_rhizo$PC2, rotated_S_rhizo$PC3)

selection_S_burkho_34 <- farthest_distance(rotated_S_burkho$PC3, rotated_S_burkho$PC4)
selection_S_para_34 <- farthest_distance(rotated_S_para$PC3, rotated_S_para$PC4)
selection_S_agro_34 <- farthest_distance(rotated_S_agro$PC3, rotated_S_agro$PC4)
selection_S_bacillus_34 <- farthest_distance(rotated_S_bacillus$PC3, rotated_S_bacillus$PC4)
selection_S_rhizo_34 <- farthest_distance(rotated_S_rhizo$PC3, rotated_S_rhizo$PC4)

selection_S_burkho_45 <- farthest_distance(rotated_S_burkho$PC4, rotated_S_burkho$PC5)
selection_S_para_45 <- farthest_distance(rotated_S_para$PC4, rotated_S_para$PC5)
selection_S_agro_45 <- farthest_distance(rotated_S_agro$PC4, rotated_S_agro$PC5)
selection_S_bacillus_45 <- farthest_distance(rotated_S_bacillus$PC4, rotated_S_bacillus$PC5)
selection_S_rhizo_45 <- farthest_distance(rotated_S_rhizo$PC4, rotated_S_rhizo$PC5)

selection_S_burkho_56 <- farthest_distance(rotated_S_burkho$PC5, rotated_S_burkho$PC6)
selection_S_para_56 <- farthest_distance(rotated_S_para$PC5, rotated_S_para$PC6)
selection_S_agro_56 <- farthest_distance(rotated_S_agro$PC5, rotated_S_agro$PC6)
selection_S_bacillus_56 <- farthest_distance(rotated_S_bacillus$PC5, rotated_S_bacillus$PC6)
selection_S_rhizo_56 <- farthest_distance(rotated_S_rhizo$PC5, rotated_S_rhizo$PC6)

min_S_burkho12 <- min(rotated_S_burkho$PC1)
min_S_rhizo12 <- min(rotated_S_rhizo$PC1)
min_S_agro12 <- min(rotated_S_agro$PC1)
min_S_para12 <- min(rotated_S_para$PC1)
min_S_bacillus12 <- min(rotated_S_bacillus$PC1)

min_S_burkho23 <- min(rotated_S_burkho$PC2)
min_S_rhizo23 <- min(rotated_S_rhizo$PC2)
min_S_agro23 <- min(rotated_S_agro$PC2)
min_S_para23 <- min(rotated_S_para$PC2)
min_S_bacillus23 <- min(rotated_S_bacillus$PC2)

min_S_burkho34 <- min(rotated_S_burkho$PC3)
min_S_rhizo34 <- min(rotated_S_rhizo$PC3)
min_S_agro34 <- min(rotated_S_agro$PC3)
min_S_para34 <- min(rotated_S_para$PC3)
min_S_bacillus34 <- min(rotated_S_bacillus$PC3)

min_S_burkho45 <- min(rotated_S_burkho$PC3)
min_S_rhizo45 <- min(rotated_S_rhizo$PC3)
min_S_agro45 <- min(rotated_S_agro$PC3)
min_S_para45 <- min(rotated_S_para$PC3)
min_S_bacillus45 <- min(rotated_S_bacillus$PC3)

min_S_burkho56 <- min(rotated_S_burkho$PC5)
min_S_rhizo56 <- min(rotated_S_rhizo$PC5)
min_S_agro56 <- min(rotated_S_agro$PC5)
min_S_para56 <- min(rotated_S_para$PC5)
min_S_bacillus56 <- min(rotated_S_bacillus$PC5)

max_S_burkho12 <- max(rotated_S_burkho$PC1)
max_S_rhizo12 <- max(rotated_S_rhizo$PC1)
max_S_agro12 <- max(rotated_S_agro$PC1)
max_S_para12 <- max(rotated_S_para$PC1)
max_S_bacillus12 <- max(rotated_S_bacillus$PC1)

max_S_burkho23 <- max(rotated_S_burkho$PC2)
max_S_rhizo23 <- max(rotated_S_rhizo$PC2)
max_S_agro23 <- max(rotated_S_agro$PC2)
max_S_para23 <- max(rotated_S_para$PC2)
max_S_bacillus23 <- max(rotated_S_bacillus$PC2)

max_S_burkho34 <- max(rotated_S_burkho$PC3)
max_S_rhizo34 <- max(rotated_S_rhizo$PC3)
max_S_agro34 <- max(rotated_S_agro$PC3)
max_S_para34 <- max(rotated_S_para$PC3)
max_S_bacillus34 <- max(rotated_S_bacillus$PC3)

max_S_burkho45 <- max(rotated_S_burkho$PC3)
max_S_rhizo45 <- max(rotated_S_rhizo$PC3)
max_S_agro45 <- max(rotated_S_agro$PC3)
max_S_para45 <- max(rotated_S_para$PC3)
max_S_bacillus45 <- max(rotated_S_bacillus$PC3)

max_S_burkho56 <- max(rotated_S_burkho$PC5)
max_S_rhizo56 <- max(rotated_S_rhizo$PC5)
max_S_agro56 <- max(rotated_S_agro$PC5)
max_S_para56 <- max(rotated_S_para$PC5)
max_S_bacillus56 <- max(rotated_S_bacillus$PC5)
# Plot the graphs

par(mfrow = c(1, 1))
par(mfrow = c(1, 5))

p_S_agro_12 <- ggplot(rotated_S_agro, aes(x = PC1, y = PC2 )) + geom_point()
p_S_agro_12 <- p_S_agro_12 + ggtitle("Agrobacterium") + xlim( min_S_agro12-40, max_S_agro12) +
  geom_text(data = rotated_S_agro[selection_S_agro_12,], 
            aes(PC1, PC2, label = df_agro$X[selection_S_agro_12]),
            hjust = 1, vjust = -.2, size = 6)  +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_agro_12

p_S_burkho_12 <- ggplot(rotated_S_burkho, aes(x = PC1, y = PC2 )) + geom_point() 
p_S_burkho_12 <- p_S_burkho_12 + ggtitle("Burkholderia") + xlim( min_S_burkho12-40, max_S_burkho12) +
  geom_text(data = rotated_S_burkho[selection_S_burkho_12,], 
            aes(PC1, PC2, label = df_burkho$X[selection_S_burkho_12]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_burkho_12

p_S_para_12 <- ggplot(rotated_S_para, aes(x = PC1, y = PC2 )) + geom_point() 
p_S_para_12 <- p_S_para_12 + ggtitle("Paracoccus") + xlim( min_S_para12-40, max_S_para12) +
  geom_text(data = rotated_S_para[selection_S_para_12,], 
            aes(PC1, PC2, label = df_para$X[selection_S_para_12]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_para_12

p_S_bacillus_12 <- ggplot(rotated_S_bacillus, aes(x = PC1, y = PC2 )) + geom_point() 
p_S_bacillus_12 <- p_S_bacillus_12 + ggtitle("Bacillus") + xlim( min_S_bacillus12-40, max_S_bacillus12) +
  geom_text(data = rotated_S_bacillus[selection_S_bacillus_12,], 
            aes(PC1, PC2, label = df_bacillus$X[selection_S_bacillus_12]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_bacillus_12

p_S_rhizo_12 <- ggplot(rotated_S_rhizo, aes(x = PC1, y = PC2 )) + geom_point() 
p_S_rhizo_12 <- p_S_rhizo_12 + ggtitle("Rhizobium") + xlim( min_S_rhizo12-40, max_S_rhizo12) +
  geom_text(data = rotated_S_rhizo[selection_S_rhizo_12,], 
            aes(PC1, PC2, label = df_rhizobium$X[selection_S_rhizo_12]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_rhizo_12

grid.arrange(p_S_agro_12, p_S_burkho_12, p_S_rhizo_12, p_S_para_12, ncol = 4)


p_S_agro_23 <- ggplot(rotated_S_agro, aes(x = PC2, y = PC3 )) + geom_point()
p_S_agro_23 <- p_S_agro_23 + ggtitle("Agrobacterium") + xlim( min_S_agro23-40, max_S_agro23) +
  geom_text(data = rotated_S_agro[selection_S_agro_23,], 
            aes(PC2, PC3, label = df_agro$X[selection_S_agro_23]),
            hjust = 1, vjust = -.2, size = 6)  +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_agro_23

p_S_burkho_23 <- ggplot(rotated_S_burkho, aes(x = PC2, y = PC3 )) + geom_point() 
p_S_burkho_23 <- p_S_burkho_23 + ggtitle("Burkholderia") + xlim( min_S_burkho23-40, max_S_burkho23) +
  geom_text(data = rotated_S_burkho[selection_S_burkho_23,], 
            aes(PC2, PC3, label = df_burkho$X[selection_S_burkho_23]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_burkho_23

p_S_para_23 <- ggplot(rotated_S_para, aes(x = PC2, y = PC3 )) + geom_point() 
p_S_para_23 <- p_S_para_23 + ggtitle("Paracoccus") + xlim( min_S_para23-40, max_S_para23) +
  geom_text(data = rotated_S_para[selection_S_para_23,], 
            aes(PC2, PC3, label = df_para$X[selection_S_para_23]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_para_23

p_S_bacillus_23 <- ggplot(rotated_S_bacillus, aes(x = PC2, y = PC3 )) + geom_point() 
p_S_bacillus_23 <- p_S_bacillus_23 + ggtitle("Bacillus") + xlim( min_S_bacillus23-40, max_S_bacillus23) +
  geom_text(data = rotated_S_bacillus[selection_S_bacillus_23,], 
            aes(PC2, PC3, label = df_bacillus$X[selection_S_bacillus_23]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_bacillus_23

p_S_rhizo_23 <- ggplot(rotated_S_rhizo, aes(x = PC2, y = PC3 )) + geom_point() 
p_S_rhizo_23 <- p_S_rhizo_23 + ggtitle("Rhizobium") + xlim( min_S_rhizo23-40, max_S_rhizo23) +
  geom_text(data = rotated_S_rhizo[selection_S_rhizo_23,], 
            aes(PC2, PC3, label = df_rhizobium$X[selection_S_rhizo_23]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_rhizo_23

grid.arrange(p_S_agro_23, p_S_burkho_23, p_S_rhizo_23, p_S_para_23, ncol = 4)


p_S_agro_34 <- ggplot(rotated_S_agro, aes(x = PC3, y = PC4 )) + geom_point()
p_S_agro_34 <- p_S_agro_34 + ggtitle("Agrobacterium") + xlim( min_S_agro34-40, max_S_agro34) +
  geom_text(data = rotated_S_agro[selection_S_agro_34,], 
            aes(PC3, PC4, label = df_agro$X[selection_S_agro_34]),
            hjust = 1, vjust = -.2, size = 6)  +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_agro_34

p_S_burkho_34 <- ggplot(rotated_S_burkho, aes(x = PC3, y = PC4 )) + geom_point() 
p_S_burkho_34 <- p_S_burkho_34 + ggtitle("Burkholderia") + xlim( min_S_burkho34-40, max_S_burkho34) +
  geom_text(data = rotated_S_burkho[selection_S_burkho_34,], 
            aes(PC3, PC4, label = df_burkho$X[selection_S_burkho_34]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_burkho_34

p_S_para_34 <- ggplot(rotated_S_para, aes(x = PC3, y = PC4 )) + geom_point() 
p_S_para_34 <- p_S_para_34 + ggtitle("Paracoccus") + xlim( min_S_para34-40, max_S_para34) +
  geom_text(data = rotated_S_para[selection_S_para_34,], 
            aes(PC3, PC4, label = df_para$X[selection_S_para_34]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_para_34

p_S_bacillus_34 <- ggplot(rotated_S_bacillus, aes(x = PC3, y = PC4 )) + geom_point() 
p_S_bacillus_34 <- p_S_bacillus_34 + ggtitle("Bacillus") + xlim( min_S_bacillus34-40, max_S_bacillus34) +
  geom_text(data = rotated_S_bacillus[selection_S_bacillus_34,], 
            aes(PC3, PC4, label = df_bacillus$X[selection_S_bacillus_34]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_bacillus_34

p_S_rhizo_34 <- ggplot(rotated_S_rhizo, aes(x = PC3, y = PC4 )) + geom_point() 
p_S_rhizo_34 <- p_S_rhizo_34 + ggtitle("Rhizobium") + xlim( min_S_rhizo34-40, max_S_rhizo34) +
  geom_text(data = rotated_S_rhizo[selection_S_rhizo_34,], 
            aes(PC3, PC4, label = df_rhizobium$X[selection_S_rhizo_34]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_rhizo_34

grid.arrange(p_S_agro_34, p_S_burkho_34, p_S_rhizo_34, p_S_para_34, ncol = 4)


p_S_agro_45 <- ggplot(rotated_S_agro, aes(x = PC4, y = PC5 )) + geom_point()
p_S_agro_45 <- p_S_agro_45 + ggtitle("Agrobacterium") + xlim( min_S_agro45-40, max_S_agro45) +
  geom_text(data = rotated_S_agro[selection_S_agro_45,], 
            aes(PC4, PC5, label = df_agro$X[selection_S_agro_45]),
            hjust = 1, vjust = -.2, size = 6)  +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_agro_45

p_S_burkho_45 <- ggplot(rotated_S_burkho, aes(x = PC4, y = PC5 )) + geom_point() 
p_S_burkho_45 <- p_S_burkho_45 + ggtitle("Burkholderia") + xlim( min_S_burkho45-40, max_S_burkho45) +
  geom_text(data = rotated_S_burkho[selection_S_burkho_45,], 
            aes(PC4, PC5, label = df_burkho$X[selection_S_burkho_45]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_burkho_45

p_S_para_45 <- ggplot(rotated_S_para, aes(x = PC4, y = PC5 )) + geom_point() 
p_S_para_45 <- p_S_para_45 + ggtitle("Paracoccus") + xlim( min_S_para45-40, max_S_para45) +
  geom_text(data = rotated_S_para[selection_S_para_45,], 
            aes(PC4, PC5, label = df_para$X[selection_S_para_45]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_para_45

p_S_bacillus_45 <- ggplot(rotated_S_bacillus, aes(x = PC4, y = PC5 )) + geom_point() 
p_S_bacillus_45 <- p_S_bacillus_45 + ggtitle("Bacillus") + xlim( min_S_bacillus45-40, max_S_bacillus45) +
  geom_text(data = rotated_S_bacillus[selection_S_bacillus_45,], 
            aes(PC4, PC5, label = df_bacillus$X[selection_S_bacillus_45]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_bacillus_45

p_S_rhizo_45 <- ggplot(rotated_S_rhizo, aes(x = PC4, y = PC5 )) + geom_point() 
p_S_rhizo_45 <- p_S_rhizo_45 + ggtitle("Rhizobium") + xlim( min_S_rhizo45-40, max_S_rhizo45) +
  geom_text(data = rotated_S_rhizo[selection_S_rhizo_45,], 
            aes(PC4, PC5, label = df_rhizobium$X[selection_S_rhizo_45]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_rhizo_45

grid.arrange(p_S_agro_45, p_S_burkho_45, p_S_rhizo_45, p_S_para_45, ncol = 4)


p_S_agro_56 <- ggplot(rotated_S_agro, aes(x = PC5, y = PC6 )) + geom_point()
p_S_agro_56 <- p_S_agro_56 + ggtitle("Agrobacterium") + xlim( min_S_agro56-40, max_S_agro56) +
  geom_text(data = rotated_S_agro[selection_S_agro_56,], 
            aes(PC5, PC6, label = df_agro$X[selection_S_agro_56]),
            hjust = 1, vjust = -.2, size = 6)  +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_agro_56

p_S_burkho_56 <- ggplot(rotated_S_burkho, aes(x = PC5, y = PC6 )) + geom_point() 
p_S_burkho_56 <- p_S_burkho_56 + ggtitle("Burkholderia") + xlim( min_S_burkho56-40, max_S_burkho56) +
  geom_text(data = rotated_S_burkho[selection_S_burkho_56,], 
            aes(PC5, PC6, label = df_burkho$X[selection_S_burkho_56]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_burkho_56

p_S_para_56 <- ggplot(rotated_S_para, aes(x = PC5, y = PC6 )) + geom_point() 
p_S_para_56 <- p_S_para_56 + ggtitle("Paracoccus") + xlim( min_S_para56-40, max_S_para56) +
  geom_text(data = rotated_S_para[selection_S_para_56,], 
            aes(PC5, PC6, label = df_para$X[selection_S_para_56]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_para_56 

p_S_bacillus_56 <- ggplot(rotated_S_bacillus, aes(x = PC5, y = PC6 )) + geom_point() 
p_S_bacillus_56 <- p_S_bacillus_56 + ggtitle("Bacillus") + xlim( min_S_bacillus56-40, max_S_bacillus56) +
  geom_text(data = rotated_S_bacillus[selection_S_bacillus_56,], 
            aes(PC5, PC6, label = df_bacillus$X[selection_S_bacillus_56]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))

p_S_bacillus_56

p_S_rhizo_56 <- ggplot(rotated_S_rhizo, aes(x = PC5, y = PC6 )) + geom_point() 
p_S_rhizo_56 <- p_S_rhizo_56 + ggtitle("Rhizobium") + xlim( min_S_rhizo56-40, max_S_rhizo56) +
  geom_text(data = rotated_S_rhizo[selection_S_rhizo_56,], 
            aes(PC5, PC6, label = df_rhizobium$X[selection_S_rhizo_56]),
            hjust = 1.25, vjust = -.2, size = 6) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_S_rhizo_56

grid.arrange(p_S_agro_56, p_S_burkho_56, p_S_rhizo_56, p_S_para_56, ncol = 4)

# Reaction plots

rotated_ST_burkho <- as.data.frame(ST_burkho.pca$x)
rotated_ST_agro <- as.data.frame(ST_agro.pca$x)
rotated_ST_para <- as.data.frame(ST_para.pca$x)
rotated_ST_bacillus <- as.data.frame(ST_bacillus.pca$x)
rotated_ST_rhizo <- as.data.frame(ST_rhizo.pca$x)

selection_ST_burkho_12 <- farthest_distance(rotated_ST_burkho$PC1, rotated_ST_burkho$PC2)
selection_ST_para_12 <- farthest_distance(rotated_ST_para$PC1, rotated_ST_para$PC2)
selection_ST_agro_12 <- farthest_distance(rotated_ST_agro$PC1, rotated_ST_agro$PC2)
selection_ST_bacillus_12 <- farthest_distance(rotated_ST_bacillus$PC1, rotated_ST_bacillus$PC2)
selection_ST_rhizo_12 <- farthest_distance(rotated_ST_rhizo$PC1, rotated_ST_rhizo$PC2)

selection_ST_burkho_23 <- farthest_distance(rotated_ST_burkho$PC2, rotated_ST_burkho$PC3)
selection_ST_para_23 <- farthest_distance(rotated_ST_para$PC2, rotated_ST_para$PC3)
selection_ST_agro_23 <- farthest_distance(rotated_ST_agro$PC2, rotated_ST_agro$PC3)
selection_ST_bacillus_23 <- farthest_distance(rotated_ST_bacillus$PC2, rotated_ST_bacillus$PC3)
selection_ST_rhizo_23 <- farthest_distance(rotated_ST_rhizo$PC2, rotated_ST_rhizo$PC3)

selection_ST_burkho_34 <- farthest_distance(rotated_ST_burkho$PC3, rotated_ST_burkho$PC4)
selection_ST_para_34 <- farthest_distance(rotated_ST_para$PC3, rotated_ST_para$PC4)
selection_ST_agro_34 <- farthest_distance(rotated_ST_agro$PC3, rotated_ST_agro$PC4)
selection_ST_bacillus_34 <- farthest_distance(rotated_ST_bacillus$PC3, rotated_ST_bacillus$PC4)
selection_ST_rhizo_34 <- farthest_distance(rotated_ST_rhizo$PC3, rotated_ST_rhizo$PC4)

selection_ST_burkho_45 <- farthest_distance(rotated_ST_burkho$PC4, rotated_ST_burkho$PC5)
selection_ST_para_45 <- farthest_distance(rotated_ST_para$PC4, rotated_ST_para$PC5)
selection_ST_agro_45 <- farthest_distance(rotated_ST_agro$PC4, rotated_ST_agro$PC5)
selection_ST_bacillus_45 <- farthest_distance(rotated_ST_bacillus$PC4, rotated_ST_bacillus$PC5)
selection_ST_rhizo_45 <- farthest_distance(rotated_ST_rhizo$PC4, rotated_ST_rhizo$PC5)

selection_ST_burkho_56 <- farthest_distance(rotated_ST_burkho$PC5, rotated_ST_burkho$PC6)
selection_ST_para_56 <- farthest_distance(rotated_ST_para$PC5, rotated_ST_para$PC6)
selection_ST_agro_56 <- farthest_distance(rotated_ST_agro$PC5, rotated_ST_agro$PC6)
selection_ST_bacillus_56 <- farthest_distance(rotated_ST_bacillus$PC5, rotated_ST_bacillus$PC6)
selection_ST_rhizo_56 <- farthest_distance(rotated_ST_rhizo$PC5, rotated_ST_rhizo$PC6)



# Write csv files with best reactions
table_ST_agro_12 <-  cbind(n = 1:length(selection_ST_agro_12), Reaction = names(df_agro)[selection_ST_agro_12+ 1], NAME.EQ = df_agro.NAME.EQ[selection_ST_agro_12])
table_ST_rhizo_12 <-  cbind(n = 1:length(selection_ST_rhizo_12), Reaction = names(df_rhizobium)[selection_ST_rhizo_12+ 1], NAME.EQ = df_rhizo.NAME.EQ[selection_ST_rhizo_12])
table_ST_para_12 <-  cbind(n = 1:length(selection_ST_para_12), Reaction = names(df_para)[selection_ST_para_12+ 1], NAME.EQ = df_para.NAME.EQ[selection_ST_para_12])
table_ST_burkho_12 <-  cbind(n = 1:length(selection_ST_burkho_12), Reaction = names(df_burkho)[selection_ST_burkho_12+ 1], NAME.EQ = df_burkho.NAME.EQ[selection_ST_burkho_12])
table_ST_bacillus_12 <-  cbind(n = 1:length(selection_ST_bacillus_12), Reaction = names(df_bacillus)[selection_ST_bacillus_12+ 1], NAME.EQ = df_bacillus.NAME.EQ[selection_ST_bacillus_12])

table_ST_agro_23 <-  cbind(n = 1:length(selection_ST_agro_23), Reaction = names(df_agro)[selection_ST_agro_23+ 1], NAME.EQ = df_agro.NAME.EQ[selection_ST_agro_23])
table_ST_rhizo_23 <-  cbind(n = 1:length(selection_ST_rhizo_23), Reaction = names(df_rhizobium)[selection_ST_rhizo_23+ 1], NAME.EQ = df_rhizo.NAME.EQ[selection_ST_rhizo_23])
table_ST_para_23 <-  cbind(n = 1:length(selection_ST_para_23), Reaction = names(df_para)[selection_ST_para_23+ 1], NAME.EQ = df_para.NAME.EQ[selection_ST_para_23])
table_ST_burkho_23 <-  cbind(n = 1:length(selection_ST_burkho_23), Reaction = names(df_burkho)[selection_ST_burkho_23+ 1], NAME.EQ = df_burkho.NAME.EQ[selection_ST_burkho_23])
table_ST_bacillus_23 <-  cbind(n = 1:length(selection_ST_bacillus_23), Reaction = names(df_bacillus)[selection_ST_bacillus_23+ 1], NAME.EQ = df_bacillus.NAME.EQ[selection_ST_bacillus_23])

table_ST_agro_34 <-  cbind(n = 1:length(selection_ST_agro_34), Reaction = names(df_agro)[selection_ST_agro_34+ 1], NAME.EQ = df_agro.NAME.EQ[selection_ST_agro_34])
table_ST_rhizo_34 <-  cbind(n = 1:length(selection_ST_rhizo_34), Reaction = names(df_rhizobium)[selection_ST_rhizo_34+ 1], NAME.EQ = df_rhizo.NAME.EQ[selection_ST_rhizo_34])
table_ST_para_34 <-  cbind(n = 1:length(selection_ST_para_34), Reaction = names(df_para)[selection_ST_para_34+ 1], NAME.EQ = df_para.NAME.EQ[selection_ST_para_34])
table_ST_burkho_34 <-  cbind(n = 1:length(selection_ST_burkho_34), Reaction = names(df_burkho)[selection_ST_burkho_34+ 1], NAME.EQ = df_burkho.NAME.EQ[selection_ST_burkho_34])
table_ST_bacillus_34 <-  cbind(n = 1:length(selection_ST_bacillus_34), Reaction = names(df_bacillus)[selection_ST_bacillus_34+ 1], NAME.EQ = df_bacillus.NAME.EQ[selection_ST_bacillus_34])

table_ST_agro_45 <-  cbind(n = 1:length(selection_ST_agro_45), Reaction = names(df_agro)[selection_ST_agro_45+ 1], NAME.EQ = df_agro.NAME.EQ[selection_ST_agro_45])
table_ST_rhizo_45 <-  cbind(n = 1:length(selection_ST_rhizo_45), Reaction = names(df_rhizobium)[selection_ST_rhizo_45+ 1], NAME.EQ = df_rhizo.NAME.EQ[selection_ST_rhizo_45])
table_ST_para_45 <-  cbind(n = 1:length(selection_ST_para_45), Reaction = names(df_para)[selection_ST_para_45+ 1], NAME.EQ = df_para.NAME.EQ[selection_ST_para_45])
table_ST_burkho_45 <-  cbind(n = 1:length(selection_ST_burkho_45), Reaction = names(df_burkho)[selection_ST_burkho_45+ 1], NAME.EQ = df_burkho.NAME.EQ[selection_ST_burkho_45])
table_ST_bacillus_45 <-  cbind(n = 1:length(selection_ST_bacillus_45), Reaction = names(df_bacillus)[selection_ST_bacillus_45+ 1], NAME.EQ = df_bacillus.NAME.EQ[selection_ST_bacillus_45])

table_ST_agro_56 <-  cbind(n = 1:length(selection_ST_agro_56), Reaction = names(df_agro)[selection_ST_agro_56+ 1], NAME.EQ = df_agro.NAME.EQ[selection_ST_agro_56])
table_ST_rhizo_56 <-  cbind(n = 1:length(selection_ST_rhizo_56), Reaction = names(df_rhizobium)[selection_ST_rhizo_56+ 1], NAME.EQ = df_rhizo.NAME.EQ[selection_ST_rhizo_56])
table_ST_para_56 <-  cbind(n = 1:length(selection_ST_para_56), Reaction = names(df_para)[selection_ST_para_56+ 1], NAME.EQ = df_para.NAME.EQ[selection_ST_para_56])
table_ST_burkho_56 <-  cbind(n = 1:length(selection_ST_burkho_56), Reaction = names(df_burkho)[selection_ST_burkho_56+ 1], NAME.EQ = df_burkho.NAME.EQ[selection_ST_burkho_56])
table_ST_bacillus_56 <-  cbind(n = 1:length(selection_ST_bacillus_56), Reaction = names(df_bacillus)[selection_ST_bacillus_56+ 1], NAME.EQ = df_bacillus.NAME.EQ[selection_ST_bacillus_56])


write.csv(table_ST_agro_12, "NAME_EQUATION_PCA/table_ST_agro_12.csv", row.names = F)
write.csv(table_ST_rhizo_12, "NAME_EQUATION_PCA/table_ST_rhizo_12.csv", row.names = F)
write.csv(table_ST_para_12, "NAME_EQUATION_PCA/table_ST_para_12.csv", row.names = F)
write.csv(table_ST_bacillus_12, "NAME_EQUATION_PCA/table_ST_bacillus_12.csv", row.names = F)
write.csv(table_ST_burkho_12, "NAME_EQUATION_PCA/table_ST_burkho_12.csv", row.names = F)

write.csv(table_ST_agro_23, "NAME_EQUATION_PCA/table_ST_agro_23.csv", row.names = F)
write.csv(table_ST_rhizo_23, "NAME_EQUATION_PCA/table_ST_rhizo_23.csv", row.names = F)
write.csv(table_ST_para_23, "NAME_EQUATION_PCA/table_ST_para_23.csv", row.names = F)
write.csv(table_ST_bacillus_23, "NAME_EQUATION_PCA/table_ST_bacillus_23.csv", row.names = F)
write.csv(table_ST_burkho_23, "NAME_EQUATION_PCA/table_ST_burkho_23.csv", row.names = F)

write.csv(table_ST_agro_34, "NAME_EQUATION_PCA/table_ST_agro_34.csv", row.names = F)
write.csv(table_ST_rhizo_34, "NAME_EQUATION_PCA/table_ST_rhizo_34.csv", row.names = F)
write.csv(table_ST_para_34, "NAME_EQUATION_PCA/table_ST_para_34.csv", row.names = F)
write.csv(table_ST_bacillus_34, "NAME_EQUATION_PCA/table_ST_bacillus_34.csv", row.names = F)
write.csv(table_ST_burkho_34, "NAME_EQUATION_PCA/table_ST_burkho_34.csv", row.names = F)

write.csv(table_ST_agro_45, "NAME_EQUATION_PCA/table_ST_agro_45.csv", row.names = F)
write.csv(table_ST_rhizo_45, "NAME_EQUATION_PCA/table_ST_rhizo_45.csv", row.names = F)
write.csv(table_ST_para_45, "NAME_EQUATION_PCA/table_ST_para_45.csv", row.names = F)
write.csv(table_ST_bacillus_45, "NAME_EQUATION_PCA/table_ST_bacillus_45.csv", row.names = F)
write.csv(table_ST_burkho_45, "NAME_EQUATION_PCA/table_ST_burkho_45.csv", row.names = F)

write.csv(table_ST_agro_56, "NAME_EQUATION_PCA/table_ST_agro_56.csv", row.names = F)
write.csv(table_ST_rhizo_56, "NAME_EQUATION_PCA/table_ST_rhizo_56.csv", row.names = F)
write.csv(table_ST_para_56, "NAME_EQUATION_PCA/table_ST_para_56.csv", row.names = F)
write.csv(table_ST_bacillus_56, "NAME_EQUATION_PCA/table_ST_bacillus_56.csv", row.names = F)
write.csv(table_ST_burkho_56, "NAME_EQUATION_PCA/table_ST_burkho_56.csv", row.names = F)


min_ST_burkho12 <- min(rotated_ST_burkho$PC1)
min_ST_rhizo12 <- min(rotated_ST_rhizo$PC1)
min_ST_agro12 <- min(rotated_ST_agro$PC1)
min_ST_para12 <- min(rotated_ST_para$PC1)
min_ST_bacillus12 <- min(rotated_ST_bacillus$PC1)

min_ST_burkho23 <- min(rotated_ST_burkho$PC2)
min_ST_rhizo23 <- min(rotated_ST_rhizo$PC2)
min_ST_agro23 <- min(rotated_ST_agro$PC2)
min_ST_para23 <- min(rotated_ST_para$PC2)
min_ST_bacillus23 <- min(rotated_ST_bacillus$PC2)

min_ST_burkho34 <- min(rotated_ST_burkho$PC3)
min_ST_rhizo34 <- min(rotated_ST_rhizo$PC3)
min_ST_agro34 <- min(rotated_ST_agro$PC3)
min_ST_para34 <- min(rotated_ST_para$PC3)
min_ST_bacillus34 <- min(rotated_ST_bacillus$PC3)

min_ST_burkho45 <- min(rotated_ST_burkho$PC3)
min_ST_rhizo45 <- min(rotated_ST_rhizo$PC3)
min_ST_agro45 <- min(rotated_ST_agro$PC3)
min_ST_para45 <- min(rotated_ST_para$PC3)
min_ST_bacillus45 <- min(rotated_ST_bacillus$PC3)

min_ST_burkho56 <- min(rotated_ST_burkho$PC5)
min_ST_rhizo56 <- min(rotated_ST_rhizo$PC5)
min_ST_agro56 <- min(rotated_ST_agro$PC5)
min_ST_para56 <- min(rotated_ST_para$PC5)
min_ST_bacillus56 <- min(rotated_ST_bacillus$PC5)

max_ST_burkho12 <- max(rotated_ST_burkho$PC1)
max_ST_rhizo12 <- max(rotated_ST_rhizo$PC1)
max_ST_agro12 <- max(rotated_ST_agro$PC1)
max_ST_para12 <- max(rotated_ST_para$PC1)
max_ST_bacillus12 <- max(rotated_ST_bacillus$PC1)

max_ST_burkho23 <- max(rotated_ST_burkho$PC2)
max_ST_rhizo23 <- max(rotated_ST_rhizo$PC2)
max_ST_agro23 <- max(rotated_ST_agro$PC2)
max_ST_para23 <- max(rotated_ST_para$PC2)
max_ST_bacillus23 <- max(rotated_ST_bacillus$PC2)

max_ST_burkho34 <- max(rotated_ST_burkho$PC3)
max_ST_rhizo34 <- max(rotated_ST_rhizo$PC3)
max_ST_agro34 <- max(rotated_ST_agro$PC3)
max_ST_para34 <- max(rotated_ST_para$PC3)
max_ST_bacillus34 <- max(rotated_ST_bacillus$PC3)

max_ST_burkho45 <- max(rotated_ST_burkho$PC3)
max_ST_rhizo45 <- max(rotated_ST_rhizo$PC3)
max_ST_agro45 <- max(rotated_ST_agro$PC3)
max_ST_para45 <- max(rotated_ST_para$PC3)
max_ST_bacillus45 <- max(rotated_ST_bacillus$PC3)

max_ST_burkho56 <- max(rotated_ST_burkho$PC5)
max_ST_rhizo56 <- max(rotated_ST_rhizo$PC5)
max_ST_agro56 <- max(rotated_ST_agro$PC5)
max_ST_para56 <- max(rotated_ST_para$PC5)
max_ST_bacillus56 <- max(rotated_ST_bacillus$PC5)
# Plot the graphs

par(mfrow = c(1, 1))
par(mfrow = c(1, 5))




p_ST_rhizo_12 <- ggplot(rotated_ST_rhizo, aes(x = PC1, y = PC2 )) + geom_point() 
p_ST_rhizo_12 <- p_ST_rhizo_12 + ggtitle("rhizobacterium") + xlim( min_ST_rhizo12-40, max_ST_rhizo12) +
  geom_text(data = rotated_ST_rhizo[selection_ST_rhizo_12,], 
            aes(PC1, PC2, label = 1:length(selection_ST_rhizo_12)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_rhizo_12

grid.arrange(p_ST_rhizo_12, tableGrob(table_ST_rhizo_12))

p_ST_para_12 <- ggplot(rotated_ST_para, aes(x = PC1, y = PC2 )) + geom_point() 
p_ST_para_12 <- p_ST_para_12 + ggtitle("parabacterium") + xlim( min_ST_para12-40, max_ST_para12) +
  geom_text(data = rotated_ST_para[selection_ST_para_12,], 
            aes(PC1, PC2, label = 1:length(selection_ST_para_12)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_para_12

grid.arrange(p_ST_para_12, tableGrob(table_ST_para_12))

p_ST_bacillus_12 <- ggplot(rotated_ST_bacillus, aes(x = PC1, y = PC2 )) + geom_point() 
p_ST_bacillus_12 <- p_ST_bacillus_12 + ggtitle("bacillusbacterium") + xlim( min_ST_bacillus12-40, max_ST_bacillus12) +
  geom_text(data = rotated_ST_bacillus[selection_ST_bacillus_12,], 
            aes(PC1, PC2, label = 1:length(selection_ST_bacillus_12)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_bacillus_12

grid.arrange(p_ST_bacillus_12, tableGrob(table_ST_bacillus_12))
p_ST_burkho_12 <- ggplot(rotated_ST_burkho, aes(x = PC1, y = PC2 )) + geom_point() 
p_ST_burkho_12 <- p_ST_burkho_12 + ggtitle("burkhobacterium") + xlim( min_ST_burkho12-40, max_ST_burkho12) +
  geom_text(data = rotated_ST_burkho[selection_ST_burkho_12,], 
            aes(PC1, PC2, label = 1:length(selection_ST_burkho_12)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_burkho_12

grid.arrange(p_ST_burkho_12, tableGrob(table_ST_burkho_12))

p_ST_agro_12 <- ggplot(rotated_ST_agro, aes(x = PC1, y = PC2 )) + geom_point() 
p_ST_agro_12 <- p_ST_agro_12 + ggtitle("agrobacterium") + xlim( min_ST_agro12-40, max_ST_agro12) +
  geom_text(data = rotated_ST_agro[selection_ST_agro_12,], 
            aes(PC1, PC2, label = 1:length(selection_ST_agro_12)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_agro_12

grid.arrange(p_ST_agro_12, tableGrob(table_ST_agro_12))


set.seed(145)
p_ST_rhizo_23 <- ggplot(rotated_ST_rhizo, aes(x = PC2, y = PC3 )) + geom_point() 
p_ST_rhizo_23 <- p_ST_rhizo_23 + ggtitle("rhizobacterium") + xlim( min_ST_rhizo23-40, max_ST_rhizo23) +
  geom_text(data = rotated_ST_rhizo[selection_ST_rhizo_23,], 
            aes(PC2, PC3, label = 1:length(selection_ST_rhizo_23)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_rhizo_23

grid.arrange(p_ST_rhizo_23, tableGrob(table_ST_rhizo_23))

p_ST_para_23 <- ggplot(rotated_ST_para, aes(x = PC2, y = PC3 )) + geom_point() 
p_ST_para_23 <- p_ST_para_23 + ggtitle("parabacterium") + xlim( min_ST_para23-40, max_ST_para23) +
  geom_text(data = rotated_ST_para[selection_ST_para_23,], 
            aes(PC2, PC3, label = 1:length(selection_ST_para_23)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_para_23

grid.arrange(p_ST_para_23, tableGrob(table_ST_para_23))

p_ST_bacillus_23 <- ggplot(rotated_ST_bacillus, aes(x = PC2, y = PC3 )) + geom_point() 
p_ST_bacillus_23 <- p_ST_bacillus_23 + ggtitle("bacillusbacterium") + xlim( min_ST_bacillus23-40, max_ST_bacillus23) +
  geom_text(data = rotated_ST_bacillus[selection_ST_bacillus_23,], 
            aes(PC2, PC3, label = 1:length(selection_ST_bacillus_23)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_bacillus_23

grid.arrange(p_ST_bacillus_23, tableGrob(table_ST_bacillus_23))
p_ST_burkho_23 <- ggplot(rotated_ST_burkho, aes(x = PC2, y = PC3 )) + geom_point() 
p_ST_burkho_23 <- p_ST_burkho_23 + ggtitle("burkhobacterium") + xlim( min_ST_burkho23-40, max_ST_burkho23) +
  geom_text(data = rotated_ST_burkho[selection_ST_burkho_23,], 
            aes(PC2, PC3, label = 1:length(selection_ST_burkho_23)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_burkho_23

grid.arrange(p_ST_burkho_23, tableGrob(table_ST_burkho_23))

p_ST_agro_23 <- ggplot(rotated_ST_agro, aes(x = PC2, y = PC3 )) + geom_point() 
p_ST_agro_23 <- p_ST_agro_23 + ggtitle("agrobacterium") + xlim( min_ST_agro23-40, max_ST_agro23) +
  geom_text(data = rotated_ST_agro[selection_ST_agro_23,], 
            aes(PC2, PC3, label = 1:length(selection_ST_agro_23)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_agro_23

grid.arrange(p_ST_agro_23, tableGrob(table_ST_agro_23))


p_ST_rhizo_34 <- ggplot(rotated_ST_rhizo, aes(x = PC3, y = PC4 )) + geom_point() 
p_ST_rhizo_34 <- p_ST_rhizo_34 + ggtitle("rhizobacterium") + xlim( min_ST_rhizo34-40, max_ST_rhizo34) +
  geom_text(data = rotated_ST_rhizo[selection_ST_rhizo_34,], 
            aes(PC3, PC4, label = 1:length(selection_ST_rhizo_34)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_rhizo_34

grid.arrange(p_ST_rhizo_34, tableGrob(table_ST_rhizo_34))

p_ST_para_34 <- ggplot(rotated_ST_para, aes(x = PC3, y = PC4 )) + geom_point() 
p_ST_para_34 <- p_ST_para_34 + ggtitle("parabacterium") + xlim( min_ST_para34-40, max_ST_para34) +
  geom_text(data = rotated_ST_para[selection_ST_para_34,], 
            aes(PC3, PC4, label = 1:length(selection_ST_para_34)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_para_34

grid.arrange(p_ST_para_34, tableGrob(table_ST_para_34))

p_ST_bacillus_34 <- ggplot(rotated_ST_bacillus, aes(x = PC3, y = PC4 )) + geom_point() 
p_ST_bacillus_34 <- p_ST_bacillus_34 + ggtitle("bacillusbacterium") + xlim( min_ST_bacillus34-40, max_ST_bacillus34) +
  geom_text(data = rotated_ST_bacillus[selection_ST_bacillus_34,], 
            aes(PC3, PC4, label = 1:length(selection_ST_bacillus_34)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_bacillus_34

grid.arrange(p_ST_bacillus_34, tableGrob(table_ST_bacillus_34))
p_ST_burkho_34 <- ggplot(rotated_ST_burkho, aes(x = PC3, y = PC4 )) + geom_point() 
p_ST_burkho_34 <- p_ST_burkho_34 + ggtitle("burkhobacterium") + xlim( min_ST_burkho34-40, max_ST_burkho34) +
  geom_text(data = rotated_ST_burkho[selection_ST_burkho_34,], 
            aes(PC3, PC4, label = 1:length(selection_ST_burkho_34)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_burkho_34

grid.arrange(p_ST_burkho_34, tableGrob(table_ST_burkho_34))

p_ST_agro_34 <- ggplot(rotated_ST_agro, aes(x = PC3, y = PC4 )) + geom_point() 
p_ST_agro_34 <- p_ST_agro_34 + ggtitle("agrobacterium") + xlim( min_ST_agro34-40, max_ST_agro34) +
  geom_text(data = rotated_ST_agro[selection_ST_agro_34,], 
            aes(PC3, PC4, label = 1:length(selection_ST_agro_34)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_agro_34

grid.arrange(p_ST_agro_34, tableGrob(table_ST_agro_34))



p_ST_rhizo_45 <- ggplot(rotated_ST_rhizo, aes(x = PC4, y = PC5 )) + geom_point() 
p_ST_rhizo_45 <- p_ST_rhizo_45 + ggtitle("rhizobacterium") + xlim( min_ST_rhizo45-40, max_ST_rhizo45) +
  geom_text(data = rotated_ST_rhizo[selection_ST_rhizo_45,], 
            aes(PC4, PC5, label = 1:length(selection_ST_rhizo_45)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_rhizo_45

grid.arrange(p_ST_rhizo_45, tableGrob(table_ST_rhizo_45))

p_ST_para_45 <- ggplot(rotated_ST_para, aes(x = PC4, y = PC5 )) + geom_point() 
p_ST_para_45 <- p_ST_para_45 + ggtitle("parabacterium") + xlim( min_ST_para45-40, max_ST_para45) +
  geom_text(data = rotated_ST_para[selection_ST_para_45,], 
            aes(PC4, PC5, label = 1:length(selection_ST_para_45)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_para_45

grid.arrange(p_ST_para_45, tableGrob(table_ST_para_45))

p_ST_bacillus_45 <- ggplot(rotated_ST_bacillus, aes(x = PC4, y = PC5 )) + geom_point() 
p_ST_bacillus_45 <- p_ST_bacillus_45 + ggtitle("bacillusbacterium") + xlim( min_ST_bacillus45-40, max_ST_bacillus45) +
  geom_text(data = rotated_ST_bacillus[selection_ST_bacillus_45,], 
            aes(PC4, PC5, label = 1:length(selection_ST_bacillus_45)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_bacillus_45

grid.arrange(p_ST_bacillus_45, tableGrob(table_ST_bacillus_45))
p_ST_burkho_45 <- ggplot(rotated_ST_burkho, aes(x = PC4, y = PC5 )) + geom_point() 
p_ST_burkho_45 <- p_ST_burkho_45 + ggtitle("burkhobacterium") + xlim( min_ST_burkho45-40, max_ST_burkho45) +
  geom_text(data = rotated_ST_burkho[selection_ST_burkho_45,], 
            aes(PC4, PC5, label = 1:length(selection_ST_burkho_45)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_burkho_45

grid.arrange(p_ST_burkho_45, tableGrob(table_ST_burkho_45))

p_ST_agro_45 <- ggplot(rotated_ST_agro, aes(x = PC4, y = PC5 )) + geom_point() 
p_ST_agro_45 <- p_ST_agro_45 + ggtitle("agrobacterium") + xlim( min_ST_agro45-40, max_ST_agro45) +
  geom_text(data = rotated_ST_agro[selection_ST_agro_45,], 
            aes(PC4, PC5, label = 1:length(selection_ST_agro_45)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_agro_45

grid.arrange(p_ST_agro_45, tableGrob(table_ST_agro_45)) 



p_ST_rhizo_56 <- ggplot(rotated_ST_rhizo, aes(x = PC5, y = PC6 )) + geom_point() 
p_ST_rhizo_56 <- p_ST_rhizo_56 + ggtitle("rhizobacterium") + xlim( min_ST_rhizo56-40, max_ST_rhizo56) +
  geom_text(data = rotated_ST_rhizo[selection_ST_rhizo_56,], 
            aes(PC5, PC6, label = 1:length(selection_ST_rhizo_56)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_rhizo_56

grid.arrange(p_ST_rhizo_56, tableGrob(table_ST_rhizo_56))

p_ST_para_56 <- ggplot(rotated_ST_para, aes(x = PC5, y = PC6 )) + geom_point() 
p_ST_para_56 <- p_ST_para_56 + ggtitle("parabacterium") + xlim( min_ST_para56-40, max_ST_para56) +
  geom_text(data = rotated_ST_para[selection_ST_para_56,], 
            aes(PC5, PC6, label = 1:length(selection_ST_para_56)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_para_56

grid.arrange(p_ST_para_56, tableGrob(table_ST_para_56))

p_ST_bacillus_56 <- ggplot(rotated_ST_bacillus, aes(x = PC5, y = PC6 )) + geom_point() 
p_ST_bacillus_56 <- p_ST_bacillus_56 + ggtitle("bacillusbacterium") + xlim( min_ST_bacillus56-40, max_ST_bacillus56) +
  geom_text(data = rotated_ST_bacillus[selection_ST_bacillus_56,], 
            aes(PC5, PC6, label = 1:length(selection_ST_bacillus_56)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_bacillus_56

grid.arrange(p_ST_bacillus_56, tableGrob(table_ST_bacillus_56))
p_ST_burkho_56 <- ggplot(rotated_ST_burkho, aes(x = PC5, y = PC6 )) + geom_point() 
p_ST_burkho_56 <- p_ST_burkho_56 + ggtitle("burkhobacterium") + xlim( min_ST_burkho56-40, max_ST_burkho56) +
  geom_text(data = rotated_ST_burkho[selection_ST_burkho_56,], 
            aes(PC5, PC6, label = 1:length(selection_ST_burkho_56)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_burkho_56

grid.arrange(p_ST_burkho_56, tableGrob(table_ST_burkho_56))

p_ST_agro_56 <- ggplot(rotated_ST_agro, aes(x = PC5, y = PC6 )) + geom_point() 
p_ST_agro_56 <- p_ST_agro_56 + ggtitle("agrobacterium") + xlim( min_ST_agro56-40, max_ST_agro56) +
  geom_text(data = rotated_ST_agro[selection_ST_agro_56,], 
            aes(PC5, PC6, label = 1:length(selection_ST_agro_56)),
            hjust = 1.25, vjust = -.2, size = 10, position=position_jitter(width=4, height=4)) +
  theme(plot.title = element_text(size=20, face="bold", vjust=2))
p_ST_agro_56

grid.arrange(p_ST_agro_56, tableGrob(table_ST_agro_56))


###
# Plotting stuff nicely
###

source("analysis 10 bacteria/highqualgraphR.R") #this requires that the file is in your working directory (check with getwd())
#Sourcing only needs to be executed once, preferably at the top of your script
highqualgraphR(p_ST_agro_12,"filename",res=1200,pointsize=12) #here x is a ggplot object


