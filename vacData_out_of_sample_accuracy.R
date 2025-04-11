
library(PLIER)

data("bloodCellMarkersIRISDMAP")
data("svmMarkers")
data("canonicalPathways")
data("vacData")
allPaths = combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)
cm=intersect(rownames(vacData), rownames(allPaths))
allPaths=allPaths[cm,]
vacData=vacData[cm,]
vacData = t(vacData)
vacData = scale(vacData)

C <- allPaths 

tCCt <- crossprod(C)
s_cross_C <- svd(tCCt)
V_C <- s_cross_C$u
D_C <- diag(as.vector(sqrt(s_cross_C$d)))
D_C_inv <- diag(as.vector(1/sqrt(s_cross_C$d)))
U_C <- C %*% V_C %*% D_C_inv
P_C <- tcrossprod(U_C)
Q_C <- diag(1, p, p) - P_C

# vacData train-test
Y <- vacData
n <- nrow(Y)

set.seed(123)
train_set <- sample(1:n, as.integer(0.8*n))
Y_train <- Y[train_set,]
Y_test <-  Y[-train_set,]

# fit model
vacData_est_k <- estimate_latent_dimension(Y_train, k_max=60)
vacData_train_fit <- compute_point_estimates(Y_train, C, k=vacData_est_k$k_hat)

vacData_train_fit <- compute_point_estimates(Y_train, C, k=50)


# use half of the genes in the test set to predict other half 
p <- ncol(Y)
set.seed(123)
impute_set <- sample(1:p, as.integer(0.5*p))

vacData_oos_predictions <- predict_oos_Y(
  Y_train, Y_test, impute_set, vacData_train_fit, P_C)

sqrt(mean((Y_test[,-impute_set] - vacData_oos_predictions$Y_pred_mean)^2))

dev.new()
plot(Y_test[,-impute_set], vacData_oos_predictions$Y_pred_mean)
abline(0, 1, col='red')











