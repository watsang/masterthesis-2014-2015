
# OBJECTIVE
# ==============
#
# 1)  Calculate SVD
# 
# 2)  Calculate diffusion kernel
#
# 
#
#



# ===========================================================


setwd("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/analysis 10 bacteria")
df_rhizobium <- read.csv("data 10 bacteria/394.7.txt")
df_agro <- read.csv("data 10 bacteria/176299.3.txt")
df_burkho <- read.csv("data 10 bacteria/269483.3.txt")
df_para <- read.csv("data 10 bacteria/318586.5.txt")
df_bacillus <- read.csv("data 10 bacteria/388400.4.txt")

ModelSEED.compounds.db <- read.csv("ModelSEED-compounds-db.csv")
ModelSEED.reactions.db <- read.csv("ModelSEED-reactions-db.csv")

find_cp <- function(stoichiometric_matrix){
  compound_names <- c()
  for (compound in 1:nrow(stoichiometric_matrix))
  {
    compound_names <- c(compound_names, toString(subset(ModelSEED.compounds.db, DATABASE == toString(stoichiometric_matrix$X[compound]) )$ABBREVIATION) )
  }
  return(compound_names)
}

toString(subset(ModelSEED.compounds.db, DATABASE == toString(df_agro$X[1]))$ABBREVIATION)

find_reaction <- function(stoichiometric_matrix){
  reaction_names <- c()
  for (reaction in 1:dim(stoichiometric_matrix)[2])
  {
    reaction_names <- append(reaction_names, toString(subset(ModelSEED.reactions.db, DATABASE == toString(names(stoichiometric_matrix)[reaction]))$NAME) )
  }
  return(reaction_names)
}


# Switch compound code with real compound name
# Switch reaction code with reaction name
names(df_rhizobium) <- c(X, find_reaction(df_rhizobium))
names(df_agro) <- c(X, find_reaction(df_agro))
names(df_burkho) <- c(X, find_reaction(df_burkho))
names(df_para) <- c(X, find_reaction(df_para))
names(df_bacillus) <- c(X, find_reaction(df_bacillus))

df_rhizobium$X

df_rhizobium$X <- find_cp(df_rhizobium)
df_agro$X <- find_cp(df_agro)
df_burkho$X <- find_cp(df_burkho)
df_para$X <- find_cp(df_para)
df_bacillus$X <- find_cp(df_bacillus)

# Calculate SVD
# --------------

df_rhizobium.svd <- svd(df_rhizobium[,2:dim(df_rhizobium)[2]])
df_agro.svd <- svd(df_agro[,2:dim(df_agro)[2]])
df_burkho.svd <- svd(df_burkho[,2:dim(df_burkho)[2]])
df_para.svd <- svd(df_para[,2:dim(df_para)[2]])
df_bacillus.svd <- svd(df_bacillus[,2:dim(df_bacillus)[2]])

beta <- 1

matrix_eigenvalues <- function(eigenvalues)
{
  zero_matrix <- matrix(data = rep(0, length(eigenvalues)^2), ncol =  length(eigenvalues), nrow = length(eigenvalues))
  diag(zero_matrix) <- eigenvalues
  return(zero_matrix)
}


# Calculate Diffusion matrices K
# -------------------------------

df_rhizobium.svd$d <- matrix_eigenvalues(df_rhizobium.svd$d)
K_rhizo <- df_rhizobium.svd$u %*% df_rhizobium.svd$d %*% t(df_rhizobium.svd$v)
names(K_rhizo) <- names(df_rhizobium)
row.names(K_rhizo) <- df_rhizobium$X

df_agro.svd$d <- matrix_eigenvalues(df_agro.svd$d)
K_agro <- df_agro.svd$u %*% df_agro.svd$d %*% t(df_agro.svd$v)
names(K_agro) <- names(df_agro)
row.names(K_agro) <- df_agro$X

df_burkho.svd$d <- matrix_eigenvalues(df_burkho.svd$d)
K_burkho <- df_burkho.svd$u %*% df_burkho.svd$d %*% t(df_burkho.svd$v)
names(K_burkho) <- names(df_burkho)
row.names(K_burkho) <- df_burkho$X

df_para.svd$d <- matrix_eigenvalues(df_para.svd$d)
K_para <- df_para.svd$u %*% df_para.svd$d %*% t(df_para.svd$v)
names(K_para) <- names(df_para)
row.names(K_para) <- df_para$X

df_bacillus.svd$d <- matrix_eigenvalues(df_bacillus.svd$d)
K_bacillus <- df_bacillus.svd$u %*% df_bacillus.svd$d %*% t(df_bacillus.svd$v)
names(K_bacillus) <- names(df_bacillus)
row.names(K_bacillus) <- df_bacillus$X



