
### run this file to see how to conduct the estimation.
library(plyr)
require(methods)

source("estimator.R") # source the estimation function
source("simulator.R") # source the simulation function

### true parameter setting
set.seed(1234)
K = 3 # set number of groups
alpha = c(0.2, 0.3, 0.5) # group ratios
theta_mat = cbind(c(0.0,0.1, 0.3, 0.5, 0.7, 1, 1.5, -1),
                  c(0.2, -0.3, 0.2, 0.1, 0.9, 0.4, -0.2, -1.5),
                  c(0.5, 0.2, 0.7, 0.2, -0.2, 1.4, -0.8, 0.5)) # true theta

### Generate the data
Z = mvrnorm(n = 200, mu = rep(0,nrow(ZSigma)), Sigma = ZSigma)  # get Z variable
#W = getBlockW(200, Nblock = 5, normalize = T) # simulate W by block model
W = getPowerLawW(N = 200, alpha = 2.5) # simulate W by power-law distribution model
Ymat_l = getY(Z, theta_mat, alpha, W, sig = 1, Time = 300)  # generate Y according to GNAR model
Ymat = Ymat_l[[1]] # get Ymat
group = Ymat_l[[2]] # get group information


### Estimate the model
em_res = EM.NAR(Ymat, W, Z, K) # EM estimation
cl_res = Cluster.NAR(Ymat, W, Z, K, method = "complete") # two-step estimation


### get the estimates for alpha theta and the group information (for EM and two-step respectively)
# for EM estimation
ind = order(em_res$alpha) # re-arrange the parameter order to be increasing
em_res$ezK = em_res$ezK[,ind]
group_em = apply(em_res$ezK, 1, which.max) 

# for two-step estimation
ind = order(cl_res$alpha)
group_cl = order(ind)[cl_res$group]

### calculate the group accuracy
c(mean(group_em!=group), mean(group_cl!=group)) 

