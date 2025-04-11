
estimate_latent_dimension <- function(Y, k_max){
  svd_Y <- svd(Y)
  n <- nrow(Y)
  jics <- sapply(1:k_max, function(x) (compute_jic(Y, svd_Y, x)))
  plot(1:k_max, jics, type='l', xlab='k', ylab='jic', main='')
  print(paste('k_hat = ', which.min(jics)))
  return(list(k_hat = which.min(jics), jics=jics, svd_Y = svd_Y))
}

compute_jic <- function(Y, svd_Y, k){
  n <- nrow(Y); p <- ncol(Y) 
  minint <- min(n ,p)
  maxint <- max(n, p)
  M <- sqrt(n)*as.matrix(svd_Y$u[,1:k])
  Lambda <- 1/sqrt(n)* as.matrix(svd_Y$v[,1:k]) %*% diag(svd_Y$d[1:k], k, k)
  
  Y_hat <- tcrossprod(M, Lambda)
  sigma_sq_hat <- colMeans((Y - Y_hat)^2) # p * 1
  tausq_est <- (mean(colSums((Y_hat)^2) / sigma_sq_hat)) / (n * k);
  
  Lambda_est <- (sqrt(n) / (n + 1/tausq_est)) * svd_Y$v[,1:k] %*% diag(svd_Y$d[1:k], k, k);
  M_Lambda_est <- sqrt(n)* svd_Y$u[,1:k] %*% t(Lambda_est)
  res <- colSums((Y - M_Lambda_est)^2) 
  sigma_sq_hat <- res / n
  
  term1 <- -n*sum(log(sqrt(sigma_sq_hat))) - 0.5*sum(res/sigma_sq_hat)
  term1 <- term1 / (n*p)
  #print(k)
  jic <- -2*term1 + (k * maxint * log(minint) / (n*p))
  #print(jic)
  return(jic)
}

compute_point_estimates <- function(
    Y, C, k = NA, k_max = 50, v0 = 1, sigma_sq0 = 1){
  
  # estimate latent dimension if unknown
  if(is.na(k)){
    k_hat <-  estimate_latent_dimension(Y, k_max)
    k <- k_hat$k_hat
    svd_Y <- k_hat$svd_Y
  }
  else{
    svd_Y <- svd(Y)
  }
  
  # extract components of svd of Y
  U <- svd_Y$u[,1:k]
  D <- diag(as.vector(svd_Y$d[1:k]))
  V <- svd_Y$v[,1:k]
  U_perp <- svd_Y$u[,-c(1:k)]
  D_perp <- diag(as.vector(svd_Y$d[-c(1:k)]))
  V_perp <- svd_Y$v[,-c(1:k)]
  
  #P_C <- C %*% solve(crossprod(C)) %*% t(C)
  tCCt <- crossprod(C)
  s_cross_C <- svd(tCCt)
  V_C <- s_cross_C$u
  D_C <- diag(as.vector(sqrt(s_cross_C$d)))
  D_C_inv <- diag(as.vector(1/sqrt(s_cross_C$d)))
  U_C <- C %*% V_C %*% D_C_inv
  P_C <- tcrossprod(U_C)
  Q_C <- diag(1, p, p) - P_C
  
  Y_P <- Y %*% P_C
  Y_Q <- Y - Y_P
  V_D <- V %*% D
  P_V_D <- P_C %*% V_D
  Q_V_D <- V_D - P_V_D 
  D_Vt_perp <- D_perp %*% t(V_perp)
  D_Vt_perp_P <-  D_Vt_perp %*% P_C
  
  # prior variances hyperparms
  tau_C <- sum((P_V_D)^2 ) / (k * sum((D_Vt_perp_P )^2) )
  print(paste('tau_C  = ', tau_C))
  tau_N <- sum((Q_V_D)^2 ) / (k * sum((D_Vt_perp - D_Vt_perp_P)^2) )
  print(paste('tau_N  = ', tau_N))
  
  # point estimates
  M <- sqrt(n) * U
  Lambda_C <- (P_V_D) * sqrt(n) / (n + 1/tau_C)
  Lambda_N <- (Q_V_D) * sqrt(n) / (n + 1/tau_N)
  sigma_sq <- (v0 * sigma_sq0 + sum((Y-tcrossprod(M, Lambda_C + Lambda_N) )^2) ) /
    (v0 + n*p -2)
  
  return(list(M=M, Lambda_C=Lambda_C, Lambda_N=Lambda_N, sigma_sq=sigma_sq,
              tau_C=tau_C, tau_N=tau_N, k=k))
}


compute_B <- function(Lambda_hat, sigma_sq_hat){
  p <- nrow(Lambda_hat)
  Lambda_hat_outer <- tcrossprod(Lambda_hat)
  B <- matrix(NA, p, p)
  for(j in 1:(p-1)){
    for(l in (j+1):p){
      B[j,l] <- sqrt(1 + (Lambda_hat_outer[j,j] * Lambda_hat_outer[l,l] + Lambda_hat_outer[j,l]) /
                       (sigma_sq_hat * (Lambda_hat_outer[j,j] + Lambda_hat_outer[l,l])))
      B[l,j] <- B[j,l]
    }
    B[j,j]  <- sqrt(1 + Lambda_hat_outer[j,j] / (2*sigma_sq_hat))
  }
  B[p,p]  <- sqrt(1 + Lambda_hat_outer[p,p] / (2*sigma_sq_hat))
}


compute_posterior_samples <- function(
    Y, Lambda_C, Lambda_N, tau_C, tau_N, sigma_sq, U_C,
    v0=1, sigma_sq_0=1, n_MC=100
){
  
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(Lambda_C)
  q <- ncol(U_C)
  
  SST <- sum(Y^2)
  SSE <- (n + 1/tau_C) * sum(Lambda_C^2) + (n + 1/tau_N) * sum(Lambda_N^2)
  sigma_sq_save <- rep(NA, n_MC)
  v_n <- v0 + n*p
  gamma_n <- v0*sigma_sq_0 + SST - SSE
  sigma_sq_save <- 1/rgamma(n_MC, v_n/2, gamma_n/2)
  
  Lambda_save <- array(NA, dim=c(p, k, n_MC))
  
  if(tau_N>tau_C){
    P_C <- tcrossprod(U_C)
    Q_C <- diag(1, p, p)- P_C
    N_C <- svd(Q_C)$u[,1:(p-q)]
  }
  
  #E <- array(rnorm(p*k*n_MC), dim=c(p, k, n_MC))
  for(t in 1:n_MC){
    if(tau_N<tau_C){
      E <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_N))
      E_2 <- U_C %*% matrix(rnorm(q*k), ncol=k) * sqrt(sigma_sq_save[t] * (1/tau_N - 1/tau_C) / (n + 1/tau_C))
    } else {
      E <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_C))
      E_2 <- N_C %*% matrix(rnorm((p-q)*k), ncol=k) * sqrt(sigma_sq_save[t] * (1/tau_C - 1/tau_N) / (n + 1/tau_N))
    }
    Lambda_save[,,t] <- Lambda_C + Lambda_N + E + E_2
  }
  return(list(Lambda_samples = Lambda_save, sigma_sq_samples = sigma_sq_save))
}


compute_posterior_samples_cc <- function(
    Y, Lambda_C, Lambda_N, tau_C, tau_N, sigma_sq, P_C,
    v0=1, sigma_sq_0=1, n_MC=100
){
  
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(Lambda_C)
  q <- ncol(U_C)
  
  Q_C <- diag(1, p, p) - P_C
  
  SST <- sum(Y^2)
  SSE <- (n + 1/tau_C) * sum(Lambda_C^2) + (n + 1/tau_N) * sum(Lambda_N^2)
  sigma_sq_save <- rep(NA, n_MC)
  v_n <- v0 + n*p
  gamma_n <- v0*sigma_sq_0 + SST - SSE
  sigma_sq_save <- 1/rgamma(n_MC, v_n/2, gamma_n/2)
  
  Lambda_save <- array(NA, dim=c(p, k, n_MC))
  
  B_C <- compute_B(Lambda_C, sigma_sq)
  rho_C <- mean(B_C[lower.tri(B_C, diag = T)])
  print(paste('rho_C = ', rho_C))
  
  
  B_N <- compute_B(Lambda_N, sigma_sq)
  rho_N <- mean(B_N[lower.tri(B_N, diag = T)])
  print(paste('rho_N = ', rho_N))
  
  #E <- array(rnorm(p*k*n_MC), dim=c(p, k, n_MC))
  for(t in 1:n_MC){
    E_C <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_C)) * rho_C
    E_C <- P_C %*% E_C
    E_N <- matrix(rnorm(p*k), ncol=k) * sqrt(sigma_sq_save[t] / (n + 1/tau_N)) * rho_N
    E_N <- Q_C %*% E_N
    Lambda_save[,,t] <- Lambda_C + Lambda_N + E_C + E_N
  }
  return(list(Lambda_samples = Lambda_save, sigma_sq_samples = sigma_sq_save))
}


samples_etas <- function(Y, Lambda_samples, sigma_sq_samples){
  
  n <- nrow(Y)
  p <- ncol(Y)
  k <- dim(Lambda_samples)[2]
  n_MC <- dim(Lambda_samples)[3]
  M_mean <- matrix(0, n, k)
  M_save <- array(0, dim=c(n, k, n_MC))
  
  for(t in 1:n_MC){
    sigma_sq <- sigma_sq_samples[t]
    Lambda <- Lambda_samples[,,t]
    prec_mat <- diag(1, k, k) + 1/sigma_sq * (crossprod(Lambda))
    svd_prec <- svd(prec_mat)
    var_sq <-  svd_prec$u %*% diag(sqrt(1/svd_prec$d)) %*% t(svd_prec$u)
    var_mat <-  svd_prec$u %*% diag(1/svd_prec$d) %*% t(svd_prec$u)
    mean_t <- 1/sigma_sq * Y %*% Lambda %*% var_mat
    M_mean <- M_mean + mean_t
    M_save[,,t] <- mean_t + matrix(rnorm(n*k), ncol=k) %*% var_sq
  }
  M_mean <- M_mean/n_MC
  
  return(list(M_mean = M_mean, M_samples = M_save))
}




compute_point_estimates_2 <- function(
    Y, C, k = NA, k_max = 50, v0 = 1, sigma_sq0 = 1,
    tol = 0.0001, iter_max = 100){
  
  # estimate latent dimension if unknown
  if(is.na(k)){
    k_hat <-  estimate_latent_dimension(Y, k_max)
    k <- k_hat$k_hat
    svd_Y <- k_hat$svd_Y
  }
  else{
    svd_Y <- svd(Y)
  }
  
  # extract components of svd of Y
  U <- svd_Y$u[,1:k]
  D <- diag(as.vector(svd_Y$d[1:k]))
  V <- svd_Y$v[,1:k]
  U_perp <- svd_Y$u[,-c(1:k)]
  D_perp <- diag(as.vector(svd_Y$d[-c(1:k)]))
  V_perp <- svd_Y$v[,-c(1:k)]
  
  
  P_C <- C %*% solve(crossprod(C)) %*% t(C)
  Q_C <- diag(1, p, p) - P_C
  
  Y_P <- Y %*% P_C
  Y_Q <- Y - Y_P
  
  V_D <- V %*% D
  P_V_D <- P_C %*% V_D
  Q_V_D <- V_D - P_V_D 
  D_Vt_perp <- D_perp %*% t(V_perp)
  D_Vt_perp_P <-  D_Vt_perp %*% P_C
  tau_C <- sum((P_V_D)^2 ) / (k * sum((D_Vt_perp_P )^2) )
  print(paste('tau_C = ', tau_C))
  tau_N <- sum((Q_V_D)^2 ) / (k * sum((D_Vt_perp - D_Vt_perp_P)^2) )
  print(paste('tau_N = ', tau_N))
  
  # initial values 
  M<- sqrt(n) * U
  Lambda_C <- (P_V_D) * sqrt(n) / (n + 1/tau_C)
  Lambda_N <- (Q_V_D) * sqrt(n) / (n + 1/tau_N)
  sigma_sq_new <- (v0 * sigma_sq0 + sum((Y-tcrossprod(M, Lambda_C + Lambda_N) )^2) ) /
    (v0 + n*p -2)
  
  it <- 1
  while(it < N_iter_max){
    print(paste0('iter : ', it))
    Lambda_C_new <- t(Y_P) %*% M %*% solve(crossprod(M) + 1/tau_C * diag(1, k, k)) 
    Lambda_N_new <- t(Y_Q) %*% M %*% solve(crossprod(M) + 1/tau_N * diag(1, k, k)) 
    var_mat <- solve(sigma_sq * diag(1, k, k) + crossprod(Lambda_C_new) + crossprod(Lambda_N_new))
    M_new <- (Y %*% (Lambda_C_new + Lambda_N_new)) %*% var_mat
    sigma_sq_new <- (
      v0 * sigma_sq0 + sum((Y-tcrossprod(M_new, Lambda_C_new + Lambda_N_new) )^2) ) /  (v0 + n*p -2)
    err <- sum((Lambda_C_new - Lambda_C)^2) + sum((Lambda_N_new - Lambda_N)^2) + sum((M_new - M)^2)
    par_old <- sum(Lambda_C^2) + sum(Lambda_N^2) + sum(M^2)
    
    print(err/par_old)
    if(err/par_old < tol){
      Lambda_C <- Lambda_C_new
      Lambda_N <- Lambda_N_new
      M <- M_new
      sigma_sq <- sigma_sq_new
      it <- N_iter_max + 1
    }
    Lambda_C <- Lambda_C_new
    Lambda_N <- Lambda_N_new
    M <- M_new
    sigma_sq <- sigma_sq_new
    it <- it + 1
  }
  
  return(list(M=M, Lambda_C=Lambda_C, Lambda_N=Lambda_N, sigma_sq=sigma_sq,
              tau_C=tau_C, tau_N=tau_N, k=k))
  
}


sample_outer_products <- function(Lambda_samples){
  p <- dim(Lambda_samples)[1]
  k <- dim(Lambda_samples)[2]
  n_MC <- dim(Lambda_samples)[3]
  
  Lambda_outer_samples <- array(NA, dim=c(p, p, n_MC))
  
  for(t in 1:n_MC){
    Lambda_outer_samples[,,t] <- tcrossprod(Lambda_samples[,,t])
  }
  return(Lambda_outer_samples)
}


