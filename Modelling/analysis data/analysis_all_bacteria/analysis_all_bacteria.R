setwd("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining")

library(Rcpp)
library(RcppEigen)
library(inline)

names_bacteria <- read.csv("analysis_all_bacteria/names_bacteria.csv")
names_bacteria
# setwd("data_organism_reactions/Stoichiometric Matrices")
names_bacteria$loc[1]
file_names <- paste("data_organism_reactions/Stoichiometric Matrices", names_bacteria$loc, sep = "/")
file_names
a <- read.table(file_names[1], row.names = 1, sep = ",", header = T)
b <- read.table(file_names[2], row.names = 1, sep = ",", header = T)

(a_new <- as.matrix(a))
str(a_new)

transCpp <- function(A)
{
  
  
}


RcppSamp <- function(X) {
  stopifnot(is.numeric(X <- as.matrix(X)),
            (nc <- ncol(X)) > 1L,
            all(X >= 0))
  .Call(transCpp, X)
}

.Call(transCpp, a_new)

RcppSamp <- function(X) {
  stopifnot(is.numeric(X <- as.matrix(X)),
            (nc <- ncol(X)) > 1L,
            all(X >= 0))
  .Call(CppSamp, X)
}

transCpp(a_new)

ftrans <- cxxfunction(signature(AA = "matrix"), transCpp, plugin = "RcppEigen")




test <- create_frames(a,b)
dim(test$S2)
check <- test[1][1]

check

create_frames <- function(S1, S2)
{
  # @param S1: stoichiometric matrix S1
  # @param S2: stoichiometric matrix S2
  # @return: returns two matrices where the compounds in both matrices are set equal to each other
  
  # Find compounds to add for each Stoichiometric matrix
  cp_add_1 <- row.names(S2)[!(row.names(S2) %in% row.names(S1))]
  cp_add_2 <- row.names(S1)[!(row.names(S1) %in% row.names(S2))]
  
  # Create zero matrix to add
  append1 <- data.frame(matrix(data=rep(0, length(cp_add_1)*ncol(S1)), nrow=length(cp_add_1), ncol=ncol(S1)), row.names=cp_add_1) 
  append2 <- data.frame(matrix(data=rep(0, length(cp_add_2)*ncol(S2)), nrow=length(cp_add_2), ncol=ncol(S2)), row.names=cp_add_2)
  
  # Add the matching reaction names for each zero matrix
  names(append1) <- names(S1)
  names(append2) <- names(S2)
  
  # Bind the zero matrices to the stoichiometric matrices
  S1 <- rbind(S1, append1)
  S2 <- rbind(S2, append2)
  
  # Order the stoichiometric matrices according to the row.names
  S1 <- S1[order(row.names(S1)),]
  S2 <- S2[order(row.names(S2)),]
  
  created_frames <- list(S1=S1, S2=S2)
  
  return(created_frames)
}

CCA <- function(S1, S2)
{
  # @param S1: stoichiometric matrix S1
  # @param S2: stoichiometric matrix S2
  # @return: returns the calculated trace lambda between S1 and S2
  
  S1 <- as.matrix(S1)
  S2 <- as.matrix(S2)
  
  S1_transposed <- t(S1)
  S2_transposed <- t(S2)
  
  term1 <- S1_transposed %*% S1
  diag(term1) <- diag(term1) + .001
  term1 <- solve(term1)
  
  term2 <- S1_transposed %*% S2
  
  term3 <- S2_transposed %*% S2
  diag(term3) <- diag(term3) + .001
  term3 <- solve(term3)
  
  term4 <- t(term2)
  
  lambda_square <- term1 %*% term2 %*% term3 %*% term4
  
  return(sum(diag(lambda_square)))
}


start <- system.time(CCA(test$S1, test$S2))

lambda <- CCA(test$S1, test$S2)
lambda
