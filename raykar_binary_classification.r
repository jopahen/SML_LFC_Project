##Implementation of Raykar-Algorithm for Binary Classification

#load auxiliary functions
source("raykar_binary_classification_functions.r")

#perform logistic regressions to get initial weights
source("logistic_regression_binary_classification.r")
print("Logistic Regressions fitted! Fitted models in model_ground_truth and model_majority.")

#load data, we work with R = 3 simulated erroneous annotators
data <- read.csv("Datasets/BreastCancerWisconsinAnnotated03.csv")[,-1]
data <- data[inds,]
N <- length(data$Diagnosis)
R <- 3

#get matrix of annotators (Ann) & predictors (X)
Ann <- as.matrix(t(cbind(data$Diagnosis1,data$Diagnosis2,data$Diagnosis3)))
X <- data[, -c(1,11,12,13,14)] #remove response variables
X <- as.matrix(X)
colnames(X) <- NULL

#initialize mu with majority vote (as commented on page 1303)
mu <- colSums(Ann) / R

#initialize sensitivities (alpha), specificities (beta) for annotators, weight vector (w) as on page 1303
alpha <- as.vector(Ann %*% mu / sum(mu))
beta <- as.vector((1 - Ann) %*% (1 - mu) / sum(1 - mu))
w <- NR_algorithm(w_init = numeric(length = length(X[1,])), mu, X, gradient, hessian)$root
theta <- c(alpha, beta, w)

#set tolerance, max_iter and counter
tol <- 10^-6
max_iter <- 100
iter <- 1

#EM-algorithm loop
print("Starting EM-algorithm to fit Raykar model...")
print("Format: number of Newton-Raphson iterations, L2-distance of parameter vector")
mu_new <- as.vector(a(alpha, Ann) * p(w, X) / (a(alpha, Ann) * p(w, X) + b(beta, Ann) *  (1 - p(w, X))))
alpha_new <- as.vector(Ann %*% mu_new / sum(mu_new))
beta_new <- as.vector((1 - Ann) %*% (1 - mu_new) / sum(1 - mu_new))
w_new <- NR_algorithm(w_init = w, mu_new, X, gradient, hessian)$root
theta_new <- c(alpha_new, beta_new, w_new)

while(sqrt(sum((theta_new - theta)^2)) > tol & iter <= max_iter){
  alpha <- alpha_new
  beta <- beta_new
  w <- w_new
  theta <- theta_new
  mu_new <- as.vector(a(alpha, Ann) * p(w, X) / (a(alpha, Ann) * p(w, X) + b(beta, Ann) *  (1 - p(w, X))))
  alpha_new <- as.vector(Ann %*% mu_new / sum(mu_new))
  beta_new <- as.vector((1 - Ann) %*% (1 - mu_new) / sum(1 - mu_new))
  NR <- NR_algorithm(w_init = w, mu_new, X, gradient, hessian)
  w_new <- NR$root
  print(NR$iterations)
  theta_new <- c(alpha_new, beta_new, w_new)
  iter <- iter + 1
  print(sqrt(sum((theta_new - theta)^2)))
}
print("Raykar-fitting complete! Estimates in raykar_out.")
raykar_out <- list(alpha = alpha_new, beta = beta_new, w = w_new, mu = mu_new)
#print(out)

#eyeball misclassification error
sum(data$Diagnosis != as.numeric(raykar_out$pred_ground_truth > 0.5))
