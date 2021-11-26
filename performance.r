library(MESS)
options(warn = -1)
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
ground_truth_test <- data_test$Diagnosis
AUCs <- numeric(4) #raykar-train, logistic-train, raykar-test, logistic-test
print("Performance Statistics:")

#create ROC-curves for Raykar model + get AUC + plot on training data + confusion matrix
par(pty = "s")
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
    sensi[i] <- sensitivity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
    speci[i] <- specificity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
}
plot(1-speci, sensi, pch = 20, col = 'blue', asp = 1,
     main = "ROC curve for Raykar (blue) & Logistic regression with majority voting (red) on training data",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'blue', lty = 3)
AUCs[1] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(raykar_out$mu), ground_truth), 1-sensitivity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(raykar_out$mu), ground_truth), specificity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Rayar (train):")
print(cmat)


#create ROC-curves for Logistic model + get AUC + plot on training data + confusion matrix
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
}
points(1-speci, sensi, col = 'red', pch = 4)
lines(1-speci, sensi, col = 'red', lty = 3)
AUCs[2] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(model_majority$fitted.values), ground_truth), 1-sensitivity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(model_majority$fitted.values), ground_truth), specificity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Logistic (train):")
print(cmat)

#create ROC-curves for Logistic model + get AUC + plot on training data + confusion matrix
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
  sensi[i] <- sensitivity(classifier(test_fitted, thresh = t[i]), ground_truth_test)
  speci[i] <- specificity(classifier(test_fitted, thresh = t[i]), ground_truth_test)
}
par(pty = "s")
plot(1-speci, sensi, col = 'red', pch = 4, asp = 1, main = "ROC curve for Raykar (blue) & Logistic regression with majority voting (red) on test data",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'red', lty = 3)
AUCs[4] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(test_fitted), ground_truth_test), 1-sensitivity(classifier(test_fitted), ground_truth_test)) * sum(ground_truth_test == 1)
neg <- c(1-specificity(classifier(test_fitted), ground_truth_test), specificity(classifier(test_fitted), ground_truth_test)) * sum(ground_truth_test == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Logistic (test):")
print(cmat)

#create ROC-curves for Raykar model + get AUC + plot on training data + confusion matrix
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
c <- sum(log((1-raykar_out$alpha)/raykar_out$beta))
logits <- log(raykar_out$alpha/(1-raykar_out$alpha)) + log(raykar_out$beta/(1-raykar_out$beta))
logits <- as.vector(logits %*% Ann_test)
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(sigmoid(Z %*% raykar_out$w + logits + c), thresh = t[i]), ground_truth_test)
  speci[i] <- specificity(classifier(sigmoid(Z %*% raykar_out$w + logits + c), thresh = t[i]), ground_truth_test)
}
points(1-speci, sensi, col = 'blue', pch = 4)
lines(1-speci, sensi, col = 'blue', lty = 3)
AUCs[3] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(sigmoid(Z %*% raykar_out$w + logits + c)), ground_truth_test), 1-sensitivity(classifier(sigmoid(Z %*% raykar_out$w + logits + c)), ground_truth_test)) * sum(ground_truth_test == 1)
neg <- c(1-specificity(classifier(sigmoid(Z %*% raykar_out$w + logits + c)), ground_truth_test), specificity(classifier(sigmoid(Z %*% raykar_out$w + logits + c)), ground_truth_test)) * sum(ground_truth_test == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Raykar (test):")
print(cmat)

options(warn = 0)

print("AUCs for Raykar (train), Logistic (train), Raykar (test), Logistic (test):")
print(AUCs)
print("ROC-plots generated!")
print("======================")
print("Complete! Warnings can be ignored...")