##Implementation of Raykar-Algorithm for Binary Classification

#load auxiliary functions
source("raykar_binary_classification_functions.r")

#perform logistic regressions to get initial weights
source("logistic_regression_binary_classification.r")
print("Logistic Regressions fitted! Fitted models in model_ground_truth and model_majority.")
print("======================")

#load data, we work with R = 3 simulated erroneous annotators
data <- read.csv("Datasets/BreastCancerWisconsinAnnotated03.csv")[,-1]
data_val <- data[-inds,]
data <- data[inds,]
N <- length(data$Diagnosis)
R <- 3

#regularization (lambda >= 0, = 0 implies no regularization)
lambda <- 0.5

#get matrix of annotators (Ann) & predictors (X)
Ann <- as.matrix(t(cbind(data$Diagnosis1,data$Diagnosis2,data$Diagnosis3)))
X <- data[, -c(1,11,12,13,14)] #remove response variables
X <- as.matrix(X)
X_val <- data_val[-c(1,11,12,13,14)]
X_val <- as.matrix(X_val)
colnames(X) <- NULL
colnames(X_val) <- NULL

#initialize mu with majority vote (as commented on page 1303)
mu_new <- colSums(Ann) / R

#initialize sensitivities (alpha), specificities (beta) for annotators, weight vector (w) as on page 1303
theta <- numeric(length = 15)
alpha_new <- as.vector(Ann %*% mu_new / sum(mu_new))
beta_new <- as.vector((1 - Ann) %*% (1 - mu_new) / sum(1 - mu_new))
w_new <- NR_algorithm(w_init = numeric(length = length(X[1,])), mu_new, X, gradient, hessian, lambda)$root
theta_new <- c(alpha_new, beta_new, w_new)
test_err <- misclas_err(as.numeric(sigmoid(X_val %*% w_new) > 0.5), data_val$Diagnosis)

#set tolerance, max_iter and counter, set validation monitor on (=TRUE) or off (=FALSE)
tol <- 10^-6
max_iter <- 100
iter <- 1
validation_monitor <- FALSE

#EM-algorithm loop
print("Starting EM-algorithm to fit Raykar model...")
print("Format: number of Newton-Raphson iterations, L2-distance of parameter vector, misclassification error on test set:")

while(sqrt(sum((theta_new - theta)^2)) > tol & iter <= max_iter){
  mu <- mu_new
  alpha <- alpha_new
  beta <- beta_new
  w <- w_new
  theta <- theta_new
  mu_new <- as.vector(a(alpha, Ann) * p(w, X) / (a(alpha, Ann) * p(w, X) + b(beta, Ann) *  (1 - p(w, X))))
  alpha_new <- as.vector(Ann %*% mu_new / sum(mu_new))
  beta_new <- as.vector((1 - Ann) %*% (1 - mu_new) / sum(1 - mu_new))
  NR <- NR_algorithm(w_init = w, mu_new, X, gradient, hessian, lambda)
  w_new <- NR$root
  print(NR$iterations)
  theta_new <- c(alpha_new, beta_new, w_new)
  iter <- iter + 1
  err <- misclas_err(as.numeric(sigmoid(X_val %*% w_new) > 0.5), data_val$Diagnosis)
  test_err <- c(test_err, err)
  print(err)
  print(sqrt(sum((theta_new - theta)^2)))
  #stop if test error stops decreasing
  if(iter > 1 & validation_monitor){
    if(test_err[iter] > test_err[iter - 1]){
      alpha_new <- alpha
      beta_new <- beta
      w_new <- w
      mu_new <- mu
      test_err <- test_err[-iter]
      break
    }
  }
}
print("Raykar-fitting complete! Estimates in raykar_out.")
raykar_out <- list(alpha = alpha_new, beta = beta_new, w = w_new, mu = mu_new,
                   fits = sigmoid(X %*% w_new), 
                   fits_val = sigmoid(X_val %*% w_new),test_err = test_err)
print("======================")