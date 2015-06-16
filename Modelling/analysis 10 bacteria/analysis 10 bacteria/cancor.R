setwd("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/analysis 10 bacteria")

library(CCA)
df_rhizobium <- read.csv("data 10 bacteria/394.7.txt")
df_agro <- read.csv("data 10 bacteria/176299.3.txt")
df_burkho <- read.csv("data 10 bacteria/269483.3.txt")
df_para <- read.csv("data 10 bacteria/318586.5.txt")
df_bacillus <- read.csv("data 10 bacteria/388400.4.txt")

# Change Lipid A disaccharide to Lipid A2
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

find_reaction <- function(stoichiometric_matrix){
  reaction_names <- c()
  for (reaction in 1:dim(stoichiometric_matrix)[2])
  {
    reaction_names <- c(reaction_names, toString(subset(ModelSEED.reactions.db, DATABASE == toString(names(stoichiometric_matrix)[reaction]))$NAME) )
  }
  return(reaction_names)
}

create_frames <- function(S1, S2)
{
  
  #   @param S1: accepts stoichiometric matrix 1
  #   @param S2: accepts stoichiometric matrix 2
  #   @return: returns two stoichiometric matrices where the rows:
  #     * consist of the same compounds
  #     * are sorted according to the index
  
  
  cp_to_add1 <- row.names(S2)[!is.element(row.names(S2), row.names(S1))]
  cp_to_add2 <- row.names(S1)[!is.element(row.names(S1), row.names(S2))]
  
  
  nadd1 <- length(cp_to_add1)
  nadd2 <- length(cp_to_add2)
  
  append1 <- matrix(rep(0, nadd1 * dim(S1)[2]), nrow = nadd1, ncol  = dim(S1)[2])
  append2 <- matrix(rep(0, nadd2 * dim(S2)[2]), nrow = nadd2, ncol  = dim(S2)[2])
  
  append1 <- as.data.frame(append1)
  append2 <- as.data.frame(append2)
  
  row.names(append1) <- cp_to_add1
  row.names(append2) <- cp_to_add2
  
  names(append1) <- names(S1)
  names(append2) <- names(S2)
  
  frame1 <- rbind(S1, append1)
  frame2 <- rbind(S2, append2)
  
  matrices <- list("S1" = frame1, "S2" = frame2)
  
  return(matrices)
  
}

dim(S_rhizo)
dim(S_burkho)

a <- create_frames(S_rhizo, S_burkho)

mdl <- estim.regul(a$S1, a$S2)

S_rhizo <- a[1]
S_burkho <- a[2]


length(names(df_rhizobium))
a <- find_reaction(df_rhizobium)
a[1]
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
S_rhizo <- data.frame( df_rhizobium[,2:dim(df_rhizobium)[2]] )
row.names(S_rhizo) <- df_rhizobium$X
names(S_rhizo) <- names(df_rhizobium)[2:dim(df_rhizobium)[2]]
dim(S_rhizo)

S_agro <- data.frame( df_agro[,2:dim(df_agro)[2]] )
row.names(S_agro) <- df_agro$X
names(S_agro) <- names(df_agro)[2:dim(df_agro)[2]]

S_burkho <- data.frame( df_burkho[,2:dim(df_burkho)[2]] )
row.names(S_burkho) <- df_burkho$X
names(S_burkho) <- names(df_burkho)[2:dim(df_burkho)[2]]

S_para <- data.frame( df_para[,2:dim(df_para)[2]] )
row.names(S_para) <- df_para$X
names(S_para) <- names(df_para)[2:dim(df_para)[2]]

S_bacillus <- data.frame( df_bacillus[,2:dim(df_bacillus)[2]] )
row.names(S_bacillus) <- df_bacillus$X
names(S_bacillus) <- names(df_bacillus)[2:dim(df_bacillus)[2]]


ST_bacillus <- t(S_bacillus)
ST_rhizo <- t(S_rhizo)
ST_agro <- t(S_agro)
ST_para <- t(S_para)
ST_burkho <- t(S_burkho)

a <- create_frames(S_rhizo, S_burkho)
a
S_rhizo <- a[1]
S_burkho <- a[2]

a[1]

dim()
dim(S_burkho)
