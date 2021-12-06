##Auxiliary functions

#logistic sigmoid
sigmoid <- function(z){
  return(1 / (1 + exp(-z)))
}

#misclassification error
misclas_err <- function(pred, label){
  return(mean(abs(pred - label)^2))
}

#functions to calculate p, a, b as on page 1302
p <- function(w, X){
  return(sigmoid(X %*% w))
}

a <- function(alpha, Ann){
  alpha_powers <- matrix(alpha ^ as.vector(Ann), ncol = 3, byrow = TRUE)
  one_minus_alpha_powers <- matrix((1-alpha) ^ as.vector(1-Ann), ncol = 3, byrow = TRUE)
  all_powers <- alpha_powers * one_minus_alpha_powers
  return(apply(all_powers,1,prod))
}

b <- function(beta, Ann){
  alpha_powers <- matrix(beta ^ as.vector(1-Ann), ncol = 3, byrow = TRUE)
  one_minus_alpha_powers <- matrix((1-beta) ^ as.vector(Ann), ncol = 3, byrow = TRUE)
  all_powers <- alpha_powers * one_minus_alpha_powers
  return(apply(all_powers,1,prod))
}

#gradient of objective (3) wrt w (page 1303)
gradient <- function(w = seq(0,1,1/8), mu, X){
  grad <- as.vector(numeric(length = length(w)))
  for(i in 1:length(mu)){
    x_i <- as.vector(X[i,])
    grad = grad + as.numeric(mu[i] - sigmoid(x_i %*% w)) * x_i
  }
  return(grad)
}

#hessian of objective (3) wrt w (page 1303)
hessian <- function(w = numeric(length = 9) + 1, X){
  hess <- matrix(numeric(length = length(w)^2), nrow = length(w))
  for(i in 1:length(X[,1])){
    x_i <- as.vector(X[i,])
    X_i <- as.matrix(X[i,] %*% t(X[i,]))
    hess = hess - as.numeric(sigmoid(x_i %*% w) * (1-sigmoid(x_i %*% w))) * X_i
  }
  return(as.matrix(hess))
}

#Newton-Raphson algorithm to approximate root of gradient of (3) wrt w (page 1303)
NR_algorithm <- function(w_init, mu, X, gradient, hessian, eta = 0.2, tol = .Machine$double.eps, maxiter = 10^3){
  w <- w_init
  w_new <- w - eta * solve(hessian(w, X)) %*% gradient(w, mu, X)
  iter <- 1
  while(sqrt(sum((w_new-w)^2)) > tol & iter <= maxiter){
    w <- w_new
    w_new <- w - eta * solve(hessian(w, X)) %*% gradient(w, mu, X)
    iter <- iter + 1
  }
  return(list(root = as.vector(w_new), gradient_value = gradient(w_new, mu, X), converged = iter <= maxiter, iterations = iter))
}