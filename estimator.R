library(MASS)

################################### Estimator Function for NAR model ##########################################

### This function is used to estimate the NAR model
betaOLS<-function(Ymat, W, Z)                                                                                  ### OLS estimation for theta: eq (2.8)
{
  Ymat1 = W%*%Ymat                                                                                             ### obtain WY
  Time = ncol(Ymat)-1                                                                                          ### Time doesn't count Y0
  if (is.null(Z))
    X = cbind(rep(1, nrow(Ymat)*Time),                                                                           ### the intercept
              as.vector(Ymat1[,-ncol(Ymat)]),                                                                    ### WY_{t-1}
              as.vector(Ymat[,-ncol(Ymat)]))
  else
    X = cbind(rep(1, nrow(Ymat)*Time),                                                                           ### the intercept
            as.vector(Ymat1[,-ncol(Ymat)]),                                                                    ### WY_{t-1}
            as.vector(Ymat[,-ncol(Ymat)]),                                                                     ### Y_{t-1}
            do.call("rbind", rep(list(Z), Time)))                                                              ### nodal covariates
  invXX = solve(crossprod(X))                                                                                  ### {t(X)X}^{-1}
  Yvec = as.vector(Ymat[,-1])                                                                                  ### the response vector
  thetaEst = invXX%*%colSums(X*Yvec)                                                                           ### estimation equation (2.8)
  sigmaHat2 = mean((Yvec - X%*%thetaEst)^2)                                                                    ### estimation for hat sigma^2
  covHat = invXX*sigmaHat2                                                                                     ### covariance for hat theta
  return(list(theta = thetaEst,
              covHat = covHat, 
              sigmaHat = sqrt(sigmaHat2)))                                                                     ### return the result
}


### EM algorithm
EM.NAR<-function(Ymat, W,  Z, K, seed = F)
{
  if (is.null(Z))
    npara = 3   # number of parameters
  else
    npara = 3+ncol(Z)
  
  alpha = rep(1/K, K)
  if (seed)
    set.seed(1234)
  
  ### set initial values for the algorithm
  nar_para = betaOLS(Ymat, W, Z)  # generate the initial estimator by the NAR model
  if (K==1)
    return(nar_para)
  theta = t(sapply(nar_para$theta, function(x) runif(K, min = x-0.05, max = x+0.05))) # generate the initial values by small perturbations of NAR estimator
  sigma = rep(nar_para$sigmaHat, K) # for simplicity set sigma the same for initial values
  
  ### get X variables
  WYmat = W%*%Ymat
  N = nrow(Ymat); Time = ncol(Ymat)
  if (is.null(Z))
    X = cbind(1, lagWY = as.vector(WYmat[,-Time]), lagY = as.vector(Ymat[,-Time]))
  else
    X = cbind(1, lagWY = as.vector(WYmat[,-Time]), lagY = as.vector(Ymat[,-Time]),
              do.call(rbind, rep(list(Z), Time - 1)))
  
  X = as.matrix(X) # all predictors
  Y = as.vector(Ymat[,-1]) # the response variable
  
  ### start the EM algorithm
  delta = 1
  ezK = matrix(0, ncol = K, nrow = nrow(Ymat))
  while(delta > 10^-5)
  {
    #cat("delta: ", delta, " ")
    
    ### E-step
    
    eps_hat = (Y - X%*%theta)%*%diag(x = 1/sigma)
    for (k in 1:K)
    {
      eps_mat = matrix(eps_hat[,k], ncol = Time - 1)
      ezK[,k] = -(Time-1)*log(sigma[k]) - rowSums(eps_mat^2/2) + log(alpha[k]) # for the stability, calculate the log-transformed probability first
    }
    
    ezK = ezK - rowMeans(ezK)
    ind = apply(ezK, 1, function(x) {
      if (any(x>500)) # if it is large enough (>500), return the largest as 1 (otherwise the R will produce an NaN)
        return(which.max(x)) # return which one is largest
      return(0)
    })
    ezK = exp(ezK)/rowSums(exp(ezK)) # exp-transform back to obtain the probability
    if (sum(ind)>0) { # code the values recorded in ind to be 0 and 1
      ii = which(ind>0)
      ezK[ii,] = 0
      ezK[cbind(ii,ind[ii])] = 1
    }
    
    ### M-step
    theta0 = theta; sigma0 = sigma; alpha0 = alpha
    for (k in 1:K)
    {
      zz = rep(ezK[,k], Time - 1)
      X_new = zz*X
      theta[,k] = ginv(crossprod(X_new, X))%*%crossprod(X_new, Y) # obtain the estimator for theta
      sigma[k] = sqrt(sum(zz*(Y - X%*%theta[,k])^2)/sum(zz)) # obtain the estimator for sigma
      
      #cat(sigma[k], " ")
      if (sigma[k]==0|is.na(sigma[k])) {
        #show(sigma[k])
        sigma[k] = 10^-4 # if the sigma is too small, change it to a small number
      }
    }
    #cat("\n")
    alpha = colSums(ezK)/N # obtain the estimator for alpha
    delta = max(abs(c(theta - theta0, sigma - sigma0, alpha - alpha0)))
  }
  return(list(theta = theta, alpha = alpha, sigma = sigma, ezK = ezK))
}


### Two step estimation

Cluster.NAR<-function(Ymat, W,  Z, K, method = "complete", plot = F, group = NULL, seed = F)
{
  N = nrow(Ymat)
  Time = ncol(Ymat)
  Ymat1 = W%*%Ymat     
  if (seed)
    set.seed(1234)
  ### obtain the regression dataset
  if (is.null(Z))
    yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                        inter = 1,
                        net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                        lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                        Y = as.vector(t(Ymat[,-1]))) 
  else
    yy_dat = data.frame(id = rep(1:N, each = Time - 1),
                        inter = 1,
                        net = as.vector(t(Ymat1[,-ncol(Ymat)])),
                        lagY = as.vector(t(Ymat[,-ncol(Ymat)])),
                        Z[rep(1:N, each = Time-1),],
                        Y = as.vector(t(Ymat[,-1])))
  
  ### for each node obtain the estimates for b_i
  paras = ddply(yy_dat, .(id), function(x){
    X = as.matrix(x[,2:4])
    invXX = ginv(crossprod(X))                                                                                   ### the response vector
    thetaEst = invXX%*%colSums(X*x$Y)   
    df = data.frame(matrix(c(thetaEst), nrow = 1))
  })
  
  colnames(paras)[-1] = c("intercept", "network", "momentum")
  
  ### scale the network and momentum parameter and calculate the nodal distances
  para_scale = apply(paras[,3:4], 2, function(x) x/max(abs(x)))#scale(paras[,3:4])
  para_dist = dist(as.matrix(para_scale))
  
  ### conduct the cluster algorithm
  if (method=="kmeans")
  {
    #nar_para = betaOLS(Ymat, W, Z)
    #ini_theta = (sapply(nar_para$theta[2:3], function(x) runif(K, min = x-0.05, max = x+0.05)))
    k_res = kmeans(para_scale, K)
    memb = k_res$cluster
  }
  else
  {
    hcl = hclust(para_dist, method = method)
    memb = cutree(hcl, k = K)
  }
  alpha = table(memb)/length(memb)
  if (plot==T)
  {
    par(mfrow = c(1,2))
    plot(paras$network, paras$momentum, col = group)
    plot(paras$network, paras$momentum, col = memb)
  }
  yy_dat$group = rep(memb, each = Time - 1) # obtain the group for each node
  ### for each group re-calculate the estimation
  theta_est = ddply(yy_dat, .(group), function(x){
    X = as.matrix(x[,2:(ncol(x)-2)])
    invXX = ginv(crossprod(X))                                                                                   ### the response vector
    thetaEst = invXX%*%colSums(X*x$Y)   
    df = data.frame(matrix(c(thetaEst), nrow = 1))
  })
  return(list(theta = t(theta_est)[-1,], alpha = alpha, group = memb))
}



### Prediction function for GNAR model (not used at this moment)


predict.GNAR<-function(Ymat, W, Z, group, theta)
{
  WYmat = W%*%Ymat
  N = nrow(Ymat); Time = ncol(Ymat)
  if (is.null(Z))
    X = cbind(1, lagWY = as.vector(WYmat[,-Time]), lagY = as.vector(Ymat[,-Time]))
  else
    X = cbind(1, lagWY = as.vector(WYmat[,-Time]), lagY = as.vector(Ymat[,-Time]),
              do.call(rbind, rep(list(Z), Time - 1)))
  Yhat = matrix((X%*%theta)[cbind(1:nrow(X), rep(group, Time-1))], nrow = N)
  R2 = 1-sum((Ymat[,-1] - Yhat)^2)/sum((Ymat[,-1] - mean(Ymat[,-1]))^2)
  p = ncol(X)-1
  adj_R2 = 1-(1-R2)*(N*(Time - 1)-1)/(N*(Time-1)-p*K-K)
  return(list(Yhat = Yhat, R2 = R2, adj_R2 = adj_R2))
}



