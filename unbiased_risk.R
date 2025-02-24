## Unbiased estimate of the risk: The following functions compute an unbiased estimate of the risk
## for an observed data Y (n by p) from N(0, Sigma0) when Sigma0 is estimated by Sigma(h) for a 
## given value of h. 

d_theta = function(lambda, n, p, h)
{
  d_vec1 = rep(0, p)
  #d_vec2 = rep(0, p)
  for(j in 1:p)
  {
    s1 = 0
    s2 = 0
    s3 = 0
    for(i in 1:p)
    {
      temp1 = (1/lambda[i] - 1/lambda[j])^2 
      temp2 =  h^2*(1/lambda[i]^2)
      s1 = s1 + (1/lambda[i])*(temp1 - temp2)/(temp1 + temp2)^2
      #s2 = s2 + h*(1/lambda[i]^2)/(temp1 + temp2)
      #s3 = s3 + h*(1/lambda[i]^2)*(1/lambda[i] - 1/lambda[j])/(temp1 + temp2)^2
    }
    d_vec1[j] = s1/p
    #d_vec2[j] = 4*(s2/p)*(s3/p)
  }
  #result = list("d_vec1" = d_vec1, "d_vec2" = d_vec2)
  return(d_vec1)
}

d_delta_inv = function(lambda, n, p, h)
{
  c1 = (1 - p/n)^2
  c2 = 2*(p/n)*(1-p/n)
  c3 = (p/n)^2
  invlambda = 1 / lambda[1:p]    # inverse of non-null eigenvalues   
  Lj = rep.row(invlambda, p)    # like 1 / lambda_j
  Lj.i = Lj - t(Lj)    # like (1 / lambda_j) - (1 / lambda_i)
  theta = rowMeans(Lj * Lj.i / (Lj.i^2 + h^2 * Lj^2))    # smoothed Stein shrinker
  #Htheta = rowMeans(Lj * (h * Lj) / (Lj.i^2 + h^2 * Lj^2)) # its conjugate
  #Atheta2 = theta^2 + Htheta^2    # its squared amplitude
  #d_mat = d_theta_Atheta(lambda, n, p, h)
  d_theta = d_theta(lambda, n, p, h)
  #d_Atheta2 = 2*theta*d_theta + d_mat$d_vec2
  d_delta_inv_vec = c1 + c2*theta + c2*(1/lambda)*d_theta #+ c3*Atheta2 + c3*(1/lambda)*d_Atheta2
  return(d_delta_inv_vec)
}


## This is the main function. It returns the estimated T1, T2 and T3. The estimate of the risk is thus
## T1 - 2*T2 + T3.
unbiased_risk = function(Y, h)
{
  n = nrow(Y)
  p = ncol(Y)
  S = (t(Y)%*%Y)/n
  eigen_S = eigen(S)
  U = eigen_S$vectors[,p:1]
  lambda = eigen_S$values[p:1]
  #qis_res = qis_delta(Y, h = h)
  res = stein_shrinker(Y, h = h)
  delta = res$deltaS
  delta_inv = 1/delta
  U = res$U
  Sig_inv = U%*%diag((1/delta))%*%t(U)
  T_1 = n*sum(lambda/(delta^2))
  A1 = (n-p-1)*Sig_inv
  A2 = sum(delta_inv)*(diag(p) - U%*%t(U))
  d_delta_inv_lambda_inv = d_delta_inv(lambda, n, p, h)
  Psi1 = delta_inv - (1/lambda)*d_delta_inv_lambda_inv
  diff_vec = rep(0, p)
  for(j in 1:p)
  {
    s = 0
    j_seq = seq(1:p)[-j]
    for(i in j_seq)
    {
      s = s + (lambda[j]/delta[j] - lambda[i]/delta[i])*(1/(lambda[j] - lambda[i]))
    }
    diff_vec[j] = s
  }
  A3 = 2*U%*%diag(Psi1)%*%t(U) + U%*%diag(diff_vec)%*%t(U)
  T_2 = sum(diag(A1 + A2 + A3))
  i_vec = rep(0, p)
  for(j in 1:p)
  {
    lm_j = lm(Y[,j] ~ Y[,-j])
    i_vec[j] = 1/summary(lm_j)$sigma^2
  }
  #T_3 = sum(i_vec*rchisq(p, n))
  T_3 = sum(n*i_vec)
  result = list("T_1" = T_1, "T_2" = T_2, "T_3" = T_3)
  return(result)
}

RSL = function(Sigma0, Sigmahat, S)
{
  Sigma0_inv = solve(Sigma0)
  Sigmahat_inv = solve(Sigmahat)
  l = sum(diag((Sigma0_inv - Sigmahat_inv)%*%(Sigma0_inv - Sigmahat_inv)%*%S))
  return(l)
}
E_RSL = function(Sigma0, R, h_seq, n)
{
  p = ncol(Sigma0)
  h0 = min((p/n)^2, 1/(p/n^2))^0.35 / p^0.35
  print(h0)
  T1_mat = matrix(0, R, length(h_seq))
  T2_mat = matrix(0, R, length(h_seq))
  T3_mat = matrix(0, R, length(h_seq))
  Sigma0_inv = solve(Sigma0)
  T4 = rep(0, length(h_seq))
  for(r in 1:R)
  {
    Y = rmvnorm(n, sigma = Sigma0)
    S = t(Y)%*%Y
    for(j in 1:length(h_seq))
    {
      #qis_res = qis_delta(Y, h = h_seq[j])
      res = stein_shrinker(Y, h = h_seq[j])
      #loss_mat[r,j] = RSL(Sigma0, qis_res$sigmahat, S)
      sigmahat = res$sigmahat
      sigmahat_inv = res$U%*%diag(1/res$deltaS)%*%t(res$U)
      #loss_mat[r,j] = sum(diag(sigmahat_inv%*%sigmahat_inv%*%S)) 
      T1_mat[r,j] = sum(diag(sigmahat_inv%*%sigmahat_inv%*%S))
      T2_mat[r,j] = sum(diag(Sigma0_inv%*%S%*%sigmahat_inv))
      T3_mat[r,j] = sum(diag(Sigma0_inv%*%Sigma0_inv%*%S))
      T4[j] = T4[j] + RSL(Sigma0, sigmahat, S)
    }
  }
  #loss_mat = T1_mat - 2*T2_mat + T3_mat
  result = list("T1" = T1_mat, "T2" = T2_mat, "T3" = T3_mat, "T4" = T4/R)
  return(result)
}


unbiased_risk_seq = function(Y, h_seq)
{
  r = rep(0, length(h_seq))
  for(i in 1:length(h_seq))
  {
    res = unbiased_risk(Y, h_seq[i])
    r[i] = res$T_1 - 2*res$T_2 + res$T_3
  }
  M = cbind(h_seq, r)
  return(M)
}