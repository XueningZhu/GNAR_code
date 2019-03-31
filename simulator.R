
### This file is for generating network structure W and responses Y.

library(MASS)
library(Matrix)
library(poweRlaw)                                                                                     ### for power-law distribution

ZSigma = 0.5^abs(outer(1:5,1:5,"-"))
######################################## useful functions #############################################

### obtain groups
get.group<-function(N, alpha)
{
  K = length(alpha) # obtain the group number
  N_group = floor(N*alpha) # calculate the number of nodes in each group
  N_group[K] = N-sum(N_group[1:(K-1)]) # re-calculate the number of nodes in the last group
  group = rep(1:K, N_group) # generate the group variable for each node
  return(group)
}

### generate the responses
getY<-function(Z, theta_mat, alpha, W, sig = 1, Time = 10)                                                 ### function for simulationg Y series given Y0, beta0, beta = (beta1, beta2), W, sigma, and T
{
  K = length(alpha) # get the number of groups
  N = nrow(W) # get the number of nodes
  group = get.group(N, alpha) # get the group for each node
  
  Ymat = matrix(0, nrow = N, ncol = Time +50 +1)                                                ### use Ymat to store the simulated data
  #Ymat[,1] = as.vector(Y0)                                                                             ### the first column of Ymat is assigned Y0
  for (i in 1:(Time+50))                 
  {
    # for each group generate the response variables
    for (k in 1:K)
    {
      beta0 = theta_mat[1,k]
      Beta = theta_mat[2:3,k]
      gamma0 = theta_mat[-(1:3),k]
      Beta0 = beta0+ Z%*%gamma0
      gg = group==k
      Ymat[gg,i+1] = as.vector(Beta0[gg] + Beta[1]*W[gg,]%*%Ymat[,i] + Beta[2]*Ymat[gg,i] 
                               + rnorm(sum(gg), sd = sig)) # generate Y according to the GNAR model
    }
  }
  return(list(Ymat = Ymat[,52:(Time+51)], group = group)) # drop the first 50 Ys (to ensure that Y arrives stationarity)
}



getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  return(as.matrix(W))
}



getPowerLawW<-function(N, alpha, normalize = T)                                                        ### get power-law network W
{
  Nfollowers = rpldis(N, 1, alpha)                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  return(W)
}
