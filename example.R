
source('helper_functions.R')


library(PLIER)
data("bloodCellMarkersIRISDMAP")
data("svmMarkers")
data("canonicalPathways")
data("vacData")

allPaths = combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)

cm=intersect(rownames(vacData), rownames(allPaths))
allPaths=allPaths[cm,]
vacData=vacData[cm,]
# we obtain C 5359x628 and Y 163x5359

vacData = t(vacData)
vacData = scale(vacData)

C <- allPaths
p <- nrow(C); q <- ncol(C); 


# compute P_C = C(C'C)^{-1}C' and Q_C = I_p - P_C
#P_C = C %*% solve(crossprod(C)) %*% t(C)
tCCt <- crossprod(C)
s_cross_C <- svd(tCCt)
V_C <- s_cross_C$u
D_C <- diag(as.vector(sqrt(s_cross_C$d)))
D_C_inv <- diag(as.vector(1/sqrt(s_cross_C$d)))
U_C <- C %*% V_C %*% D_C_inv
P_C <- tcrossprod(U_C)
Q_C <- diag(1, p, p) - P_C

n <- 300
k <- 10
sigma_sq_0 <- 2

# generate data

## loadings
Gamma_0 <- matrix(rnorm(q*k, 0, 0.1), ncol=k)
C_Gamma_0 <- C %*% Gamma_0
Psi_0 <- Q_C %*% matrix(rnorm(p*k, 0, 0.2), ncol=k)
sum(C_Gamma_0^2)/sum(Psi_0^2)
Lambda_0 <- C_Gamma_0 + Psi_0
boxplot(rowSums(Lambda_0^2) / sigma_sq_0) # snr across genes

## latent factors
M_0 <-  matrix(rnorm(n*k, 0, 1), ncol=k)

## data
Y <- M_0 %*% t(Lambda_0) +
  sqrt(sigma_sq_0) * matrix(rnorm(n*p), nrow=n) 

## true covariance components
Lambda_0_outer <- tcrossprod(Lambda_0)
Lambda_C_0_outer <- tcrossprod(C_Gamma_0)
Lambda_N_0_outer <- tcrossprod(Psi_0)

# fit point estimates
fit <- compute_point_estimates(Y, C, k=k)
Lambda_hat <- fit$Lambda_C + fit$Lambda_N

print(norm(tcrossprod(fit$Lambda_C) - Lambda_C_0_outer, type='F') / norm(Lambda_C_0_outer, type='F'))
print(norm(tcrossprod(fit$Lambda_N) - Lambda_N_0_outer, type='F') / norm(Lambda_N_0_outer, type='F'))
print(norm(tcrossprod(Lambda_hat) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F'))


# fable 
svd_Y <- svd(Y)
D <- diag(svd_Y$d[1:k])
D_perp <- diag(svd_Y$d[-c(1:k)])
V <- svd_Y$v[,1:k]
tau_L <- sum((D)^2 ) / (k * sum((D_perp)^2) )
Lambda_hat_fable <- V %*% D *  sqrt(n) / (n + 1/tau_L)
print(norm(tcrossprod(Lambda_hat_fable) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F'))


# UQ 

## sample loadings given M = M_hat
params_posterior_samples <- compute_posterior_samples_cc_2(
  Y, fit$Lambda_C, fit$Lambda_N, fit$tau_C, fit$tau_N, fit$sigma_sq, P_C
)

# latent factors
M_posterior_samples <- samples_etas(Y, params_posterior_samples$Lambda_samples, params_posterior_samples$sigma_sq_samples)

## posterior mean of latent factors
M_post_mean <- M_posterior_samples$M_mean 

## PCA estimator of latent factors
M_hat <- fit$M

## Procrustes RMSE comparison
library(vegan)

### post mean
p_mean <- procrustes(M_0, M_post_mean)
sqrt(sum((M_0 - p_mean$Yrot)^2) / (n*k))

### PCA estimator
p_hat <- procrustes(M_0, M_hat)
sqrt(sum((M_0 - p_hat$Yrot)^2) / (n*k))

# covariance components
subsample_index <- 1:500

Lambda_outer_posterior_samples <- sample_Lambda_outer(
  params_posterior_samples$Lambda_samples[subsample_index,,])

Lambda_outer_qs <- apply(Lambda_outer_posterior_samples, c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
mean((Lambda_outer_qs[1,,]<Lambda_0_outer[subsample_index, subsample_index]) &
       (Lambda_outer_qs[2,,]>Lambda_0_outer[subsample_index, subsample_index]))


# estimate k
est_k <- estimate_latent_dimension(Y, k_max=20)




