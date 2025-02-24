# EBayes-Integration
The functions in the file github_codes.R contain R functions to implement the Empricial Bayes data integration method. Specifically, let $X$ be a $N \times p$ design matrix and $y^{(t)}$ be the response vector for the $t$-th tissue where $t = 1, \ldots, n$. We assume $y^{(t)} = X\beta^{(t)} + \epsilon^{(t)}$ where $\epsilon^{(t)} \sim \mathrm{N}(0,\sigma^2)$. Suppose that $N > p$, so that ordinary least squares estimates $\hat{\beta}^{(t)}$ exist and $\hat{\beta}^{(t)} \sim \mathrm{N}(\beta^{(t)}, \sigma^2 (X'X)^{-1})$. An estimate of $\sigma^2$ can be obtained from the $t$-th tissue as $||y - X\hat{\beta}^{(t)}||_2^2/(n-p)$. Define $\hat{\sigma}^2$ as the average of these estimators and let $\Omega = \hat{\sigma}^2 (X'X)^{-1}$. Write $\hat{B}$ as the matrix with $n$ rows and $p$ columns with each row set to $\hat{\beta}^{(t)}$. Define $B$ as a matrix with $\beta^{(t)}$ as its rows. This work proposes an asymptotically optimal linear shrinkage estimator of $B$. 

The key idea is to estimate the covariance matrix of $\hat{\beta}^{(t)}$ which is then used to compute a linear shrinkage estimator of $B$ as $\widetilde{B}=B(I - C)\hat{B}$. The estimator is designed to work under nominal structural assumption on $B$: if a structure is present, then it works comparably with estimators that assume the correct structure, however, if the structural assumption is wrong, the estimator providew robust analysis.

### Typical example ###
```
library(mvtnorm)
library(pracma)
source("unbiased_risk.R")
source("github_codes.R")
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
norm(t(beta0) - Bhat2, "F") ## error for linear shrinkage estimator with default smoothing parameter 
norm(t(beta0) - Bhat1, "F") ## error for linear shrinkage estimator with smoothing parameter chosen using SURE

```
