data(Paracou618)

# Prepare the similarity matrix
DistanceMatrix <- as.matrix(Paracou618.dist)
# Similarity can be 1 minus normalized distances between species
Z <- 1 - DistanceMatrix/max(DistanceMatrix)
# Calculate diversity of order 2
Dqz(Paracou618.MC$Ns, 0, Z)

length(Paracou618.MC$Ns[Paracou618.MC$Ns > 0])

Paracou618.MC$Ns

Dqz_wk <- function (Ps, q = 1, Z = diag(length(Ps)), CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  if (is.null(colnames(Z)) | is.null(names(Ps))) {
    if (ncol(as.matrix(Z)) != length(Ps)) 
      stop("The matrix dimension must equal the probability vector length.")
    Z <- as.matrix(Z)[Ps != 0, Ps != 0]
    Ps <- Ps[Ps != 0]
  }
  else {
    Ps <- Ps[Ps != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0) 
      stop("Some species are missing in the similarity matrix.")
    Z <- as.matrix(Z)[names(Ps), names(Ps)]
  }
  Zp <- Z %*% Ps
  print(Zp)
  if (q == 1) {
    Diversity <- exp(-Ps %*% log(Zp))
  }
  else {
    Zpqm1 <- Zp^(q - 1)
    print(Zpqm1)
    Diversity <- (Ps %*% Zpqm1)^(1/(1 - q))
  }
  return(as.numeric(Diversity))
}

Dqz_wk(Ps, q = 0, Z = sim)




Ps <- as.numeric(proportions[1,])
names(Ps) <- NULL
as.numeric(Ps)

diag(sim) <- rep(1,10)
sim

proportions <- read.csv("D:/univ/2014-2015/thesis/KERMIT/bacterial similarity/data mining/proportions.csv")

sim <- matrix(rep(0,100), ncol = 10, nrow=10)
sim <- as.data.frame(sim)
names(sim) <- names(proportions)
row.names(sim) <- names(proportions)
sim


