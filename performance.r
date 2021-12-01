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
AUCs <- numeric(6) #raykar ground_truth, majority ground_truth, raykar_classifier-train, majority_classifier-train, raykar_classifier-test, majority_classifier-test
print("Performance Statistics:")

#create ROC-curves for Raykar model ground truth estimation on training data + get AUC + plot + confusion matrix
par(pty = "s")
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
    sensi[i] <- sensitivity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
    speci[i] <- specificity(classifier(raykar_out$mu, thresh = t[i]), ground_truth)
}
plot(1-speci, sensi, pch = 20, col = 'blue', asp = 1,
     main = "ROC curve for ground truth estimation on training data: \n Raykar (blue) & majority voting (red)",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'blue', lty = 3)
AUCs[1] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(raykar_out$mu), ground_truth), 1-sensitivity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(raykar_out$mu), ground_truth), specificity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Rayar (ground truth):")
print(cmat)

#create ROC-curves for majority ground truth est. on training data + get AUC + plot on training data + confusion matrix
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(Diagnosis_majority_probs, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(Diagnosis_majority_probs, thresh = t[i]), ground_truth)
}
sensi <- c(1,sensi)
speci <- c(0,speci)
points(1-speci, sensi, col = 'red', pch = 4)
lines(1-speci, sensi, col = 'red', lty = 3)
AUCs[2] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(Diagnosis_majority_probs), ground_truth), 1-sensitivity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(Diagnosis_majority_probs), ground_truth), specificity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Majority (ground truth):")
print(cmat)

#create ROC-curves for Raykar classifier on training data + get AUC + plot + confusion matrix
par(pty = "s")
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(raykar_out$fits, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(raykar_out$fits, thresh = t[i]), ground_truth)
}
plot(1-speci, sensi, pch = 20, col = 'blue', asp = 1,
     main = "ROC curve for diagnosis estimation on training data \n Raykar (blue) & majority voting (red)",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'blue', lty = 3)
AUCs[3] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(raykar_out$fits), ground_truth), 1-sensitivity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(raykar_out$fits), ground_truth), specificity(classifier(raykar_out$mu), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Rayar (classifier):")
print(cmat)

#create ROC-curves for Logistic classifier on training data + get AUC + plot + confusion matrix
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
}
points(1-speci, sensi, col = 'red', pch = 4)
lines(1-speci, sensi, col = 'red', lty = 3)
AUCs[4] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(model_majority$fitted.values), ground_truth), 1-sensitivity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 1)
neg <- c(1-specificity(classifier(model_majority$fitted.values), ground_truth), specificity(classifier(model_majority$fitted.values), ground_truth)) * sum(ground_truth == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Majority (classifier):")
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
plot(1-speci, sensi, col = 'red', pch = 4, asp = 1, main = "ROC curve for Raykar classifier (blue) & \n Majority classifier (red) on test data",
     xlab = "FPR = 1 - Specificity", ylab = "Sensitivity")
lines(1-speci, sensi, col = 'red', lty = 3)
AUCs[6] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(test_fitted), ground_truth_test), 1-sensitivity(classifier(test_fitted), ground_truth_test)) * sum(ground_truth_test == 1)
neg <- c(1-specificity(classifier(test_fitted), ground_truth_test), specificity(classifier(test_fitted), ground_truth_test)) * sum(ground_truth_test == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Majority (classifier, test):")
print(cmat)

#create ROC-curves for Raykar model + get AUC + plot on training data + confusion matrix
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
c <- sum(log((1-raykar_out$alpha)/raykar_out$beta))
logits <- log(raykar_out$alpha/(1-raykar_out$alpha)) + log(raykar_out$beta/(1-raykar_out$beta))
logits <- as.vector(logits %*% Ann_test)
logits <- Z %*% raykar_out$w # + logits + c
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(sigmoid(logits), thresh = t[i]), ground_truth_test)
  speci[i] <- specificity(classifier(sigmoid(logits), thresh = t[i]), ground_truth_test)
}
points(1-speci, sensi, col = 'blue', pch = 20)
lines(1-speci, sensi, col = 'blue', lty = 3)
AUCs[5] <- auc(1-speci, sensi)
pos <- c(sensitivity(classifier(sigmoid(logits)), ground_truth_test), 1-sensitivity(classifier(sigmoid(logits)), ground_truth_test)) * sum(ground_truth_test == 1)
neg <- c(1-specificity(classifier(sigmoid(logits)), ground_truth_test), specificity(classifier(sigmoid(logits)), ground_truth_test)) * sum(ground_truth_test == 0)
cmat <- matrix(c(pos,neg), 2)
print("Confusion Matrix Raykar (classifier, test):")
print(cmat)

#plot test-error decay during training
plot(raykar_out$test_err, col = 'orange', pch = 20, xlab = "EM-iterations", ylab = "avg. misclas. error",
     main = "Misclassification error decay \n (during training on test set)")
lines(raykar_out$test_err, col = 'orange', lty = 3)

options(warn = 0)

print("AUCs for Raykar estimate of ground truth on training data, Majority estimate of ground truth on training data, Raykar classifier (train), Majority classifier (train), Raykar classifier (test), Majority classifier (test):")
print(AUCs)
print("ROC-plots generated!")
print("======================")
print("Complete! Warnings can be ignored...")