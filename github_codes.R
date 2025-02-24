### Functions to implement EBayesIntegration ###
library(pracma)
library(mvtnorm)
rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

## The smoothed Stein Shrinker for estimating covariance matrices ##
stein_shrinker <- function(Y, k = -1, h) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  sample <- (t(Y) %*% Y) / n    # sample covariance matrix    
  spectral <- eigen(sample)    # spectral decompositon
  lambda <- spectral$values[p:1]    # sort eigenvalues in ascending order
  u <- spectral$vectors[,p:1]    # eigenvectors follow their eigenvalues
  #h <- min(c^2, 1/c^2)^0.35 / p^0.35    # smoothing parameter
  invlambda <- 1 / lambda[max(1, p-n+1):p]    # inverse of non-null eigenvalues   
  Lj <- rep.row(invlambda, min(p, n))    # like 1 / lambda_j
  Lj.i <- Lj - t(Lj)    # like (1 / lambda_j) - (1 / lambda_i)
  theta <- rowMeans(Lj * Lj.i / (Lj.i^2 + h^2 * Lj^2)) 
  theta <- sqrt(theta^2)
  # smoothed Stein shrinker
  if(p<=n)
  {
    delta <- 1/((1-c)*invlambda + 2*c*invlambda*theta) ## Optimally shrunk eigenvalues
  }else {
    delta0 <- 1 / ((c - 1) * mean(invlambda)) 
    delta <- 1/((c-1)*invlambda + 2*theta*invlambda)
    delta <- c(rep(delta0, p - n), delta)
  }
  
  
  
  deltaS <- delta * (sum(lambda) / sum(delta))    # preserve trace
  sigmahat <- u %*% diag(deltaS) %*% t(u) #reconstruct covariance matrix
  sigmahat_inv <- u%*%diag(1/deltaS)%*%t(u)
  result = list("sigmahat" = sigmahat, "sigmahat_inv" = sigmahat_inv,  "h" = h, "deltaS" = deltaS, "U" = u)
}

## Ledoit and Wolf default choice of h ##
default_h = function(n, p)
{
  c <- p / n
  h0 <- min(c^2, 1/c^2)^0.35 / p^0.35    # smoothing parameter
  return(h0)
}

## Optimal linear shirnkage of mean matrix given smoothing parameter h and covariance Omega ##
## In EBayes integration each row of X represents estimated coefficient vector from one tissue ##
## Omega is equal to sigma^2(D'D)^{-1} where D is n by p known design matrix and sigma^2 is an estimate of the noise variance ##
opt_linear_shrinkage = function(X, h, Omega)
{
  p = ncol(X)
  n = nrow(X)
  sqrt_Omega = sqrtm(Omega)
  Omega_half = sqrt_Omega$B
  Omega_half_inv = sqrt_Omega$Binv
  X_star = X%*%Omega_half_inv
  res = stein_shrinker(X_star, k = -1, h)
  Sigma_hat_inv = res$sigmahat_inv
  theta_star_hat = X_star%*%(diag(p) - Sigma_hat_inv)
  theta_hat = theta_star_hat%*%Omega_half
  return(theta_hat)
}

## Optimal linear shirnkage of mean matrix covariance Omega but smoothing parameter is chosen according to SURE ##
opt_linear_shrinkage_unb = function(X, Omega)
{
  p = ncol(X)
  n = nrow(X)
  sqrt_Omega = sqrtm(Omega)
  Omega_half = sqrt_Omega$B
  Omega_half_inv = sqrt_Omega$Binv
  X_star = X%*%Omega_half_inv
  h_seq = seq(0.01, 0.5, by = 0.02)
  h_risk = matrix(0, length(h_seq), 3)
  unbiased_risk_estimate = rep(0, length(h_seq))
  for(ii in 1:length(h_seq)) {
    res_i = unbiased_risk(Y, h_seq[ii])
    h_risk[ii, 1] = res_i$T_1
    h_risk[ii, 2] = res_i$T_2
    h_risk[ii, 3] = res_i$T_3
    unbiased_risk_estimate[ii] = h_risk[ii,1] - 2*h_risk[ii,2] + h_risk[ii,3]
  }
  
  unbiased_risk_estimate2=unbiased_risk_estimate
  unbiased_risk_estimate2[unbiased_risk_estimate<0]=Inf
  h_unb=h_seq[which.min(unbiased_risk_estimate2)]
  res = stein_shrinker(X_star, h = h_unb)
  
  Sigma_hat_inv = res$sigmahat_inv
  theta_star_hat = X_star%*%(diag(p) - Sigma_hat_inv)
  theta_hat = theta_star_hat%*%Omega_half
  return(theta_hat)
}

## Function to compute local linear shrinkage estimator ##
local_linear_shrinkage = function(X, h, K, Omega, nmcmc, burnin)
{
  n_tissues = nrow(X)
  p = ncol(X)
  Z = rmultinom(n_tissues, 1, rep(1/K, K))
  
  
  sqrt_Omega = sqrtm(Omega)
  Omega_half = sqrt_Omega$B
  Omega_half_inv = sqrt_Omega$Binv
  X_star = X%*%Omega_half_inv
  
  Sigma_inv = array(0, c(p, p, K, nmcmc))
  Sigma = array(0, c(p, p, K, nmcmc))
  B_samples1 = array(0, c(n_tissues, p, nmcmc))
  B_samples2 = array(0, c(n_tissues, p, nmcmc))
  for(ii in 1:nmcmc)
  {
    ## Update Sigma given Z ##
    for(k in 1:K)
    {
      if(sum(Z[k,])>1){
        X_k = X_star[Z[k,] == 1,]
        res_k = stein_shrinker(X_k, h = h) 
        Sigma_inv[,, k, ii] = res_k$sigmahat_inv 
        Sigma[,,k, ii] = res_k$sigmahat 
        Sigma[,,k, ii] = 0.5*(Sigma[,,k, ii] + t(Sigma[,,k, ii]))
      }else{
        Sigma[,,k, ii] = Sigma[,,k, ii-1]
        Sigma_inv[,,k, ii] = Sigma_inv[,,k, ii-1]
      }
    }
    ## Update Z given Sigma ##
    #pi_vec = c(n_tissues - sum(Z), sum(Z))/n_tissues
    pi_vec = rowSums(Z)/n_tissues
    for(i in 1:n_tissues)
    {
      probs = rep(0, K)
      for(k in 1:K)
      {
        probs[k] = pi_vec[k]*dmvnorm(X_star[i,], rep(0, p), Sigma[, , k, ii])
      }
      probs = probs/sum(probs)
      print(probs)
      for(k in 1:K)
      {
        B_samples1[i, , ii] = B_samples1[i, , ii] + probs[k]*(diag(p) - Sigma_inv[,,k,ii])%*%X_star[i,]
      }
      Z[,i] = rmultinom(1, 1, probs)    
    }
    if(ii%%100 == 0) print(ii)
    B_samples2[, , ii] = B_samples1[,,ii]%*%Omega_half
  }
  
  return(B_samples2)
}

### Example ###
N = 100
p = 20
ntissues = 40
rho_X = 0.5
Sigma_X = (1 - rho_X)*diag(p) + rho_X*matrix(1, p, p)
X = rmvnorm(N, sigma = Sigma_X)
r = 10
beta0 = matrix(rnorm(p*r), p, r)%*%matrix(rnorm(r*ntissues), r, ntissues)
sigma0 = 1
Y = matrix(0, N, ntissues)
beta = matrix(0, p, ntissues)
sigma_hat = rep(0, ntissues)
for(j in 1:ntissues)
{
  Y[,j] = X%*%beta0[,j] + sigma0*rnorm(N)
  lm_j = lm(Y[,j]~-1+X)
  beta[,j] = lm_j$coefficients
  sigma_hat[j] = summary(lm_j)$sigma
}
Omega = mean(sigma_hat)^2*solve(t(X)%*%X)
Bhat1 = opt_linear_shrinkage_unb(t(beta), Omega)
Bhat2 = opt_linear_shrinkage(t(beta), default_h(ntissues, p), Omega)
norm(t(beta0) - t(beta), "F") ## OLS error
norm(t(beta0) - Bhat1, "F") 
norm(t(beta0) - Bhat2, "F")
