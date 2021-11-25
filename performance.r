#fit Raykar and Logistic Regression models with previous routines
source("raykar_binary_classification.r")

#function to evaluate classifier based on different model fitted values and thresholds
classifier <- function(fitted_vals, thresh = 0.5){
  return(as.numeric(fitted_vals > thresh))
}

#functions to calculate specificities/sensitivities
sensitivity <- function(classified_vals, ground_truth){
  den <- length(ground_truth[ground_truth == 1])
  num <- length(classified_vals[(classified_vals == 1) & (ground_truth == 1)])
  return(num/den)
}

specificity <- function(classified_vals, ground_truth){
  den <- length(ground_truth[ground_truth == 0])
  num <- length(classified_vals[(classified_vals == 0) & (ground_truth == 0)])
  return(num/den)
}

ground_truth <- data_ground_truth$Diagnosis

#create ROC-curves for Raykar model + plot on training data
par(pty = "s")
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
    sensi[i] <- sensitivity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
    speci[i] <- specificity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
}
plot(1-speci, sensi, pch = 20, col = 'blue', asp = 1,
     main = "ROC curve for Raykar (blue) & Logistic regression with majority voting (red)",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'blue', lty = 3)

#create ROC-curves for Logistic model + plot on training data
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
}
points(1-speci, sensi, col = 'red', pch = 4)
lines(1-speci, sensi, col = 'red', lty = 3)

#create ROC-curves for Logistic model + plot on testing data
Ann_test <- as.matrix(t(cbind(data_test$Diagnosis1,data_test$Diagnosis2,data_test$Diagnosis3)))
Diagnosis_majority <- colSums(Ann_test) / R
data_test <- data_test[,-c(1, 11, 12, 13, 14)]
data_test$DiagnosisMajority <- Diagnosis_majority
Z <- as.matrix(data_test[,-10])
test_fitted <- sigmoid(Z %*% model_majority$coefficients)


t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(test_fitted, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(test_fitted, thresh = t[i]), ground_truth)
}
par(pty = "s")
plot(1-speci, sensi, col = 'red', pch = 4, asp = 1, main = "ROC curve for Raykar (blue) & Logistic regression with majority voting (red)",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'red', lty = 3)

#create ROC-curves for Raykar model + plot on testing data
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(sigmoid(Z %*% raykar_out$w), thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(sigmoid(Z %*% raykar_out$w), thresh = t[i]), ground_truth)
}
points(1-speci, sensi, col = 'blue', pch = 4)
lines(1-speci, sensi, col = 'blue', lty = 3)
