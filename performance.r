#fit Raykar and Logistic Regression models with previous routines
source("raykar_binary_classification.r")

#function to evaluate classifier based on different model fitted values and thresholds
classifier <- function(fitted_vals, thresh = 0.5){
  return(as.numeric(fitted_vals > thresh))
}

#functions to calculate specificities/sensitivities
specificity <- function(classified_vals, ground_truth){
  den <- length(ground_truth[ground_truth == 1])
  num <- length(classified_vals[(classified_vals == 1) & (ground_truth == 1)])
  return(num/den)
}

sensitivity <- function(classified_vals, ground_truth){
  den <- length(ground_truth[ground_truth == 0])
  num <- length(classified_vals[(classified_vals == 0) & (ground_truth == 0)])
  return(num/den)
}

ground_truth <- data_ground_truth$Diagnosis

#create ROC-curves for Raykar model + plot
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

#create ROC-curves for Logistic model + plot
t <- seq(0,1,0.02)
sensi <- numeric(length(t))
speci <- numeric(length(t))
for(i in 1:length(t)){
  sensi[i] <- sensitivity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
  speci[i] <- specificity(classifier(model_majority$fitted.values, thresh = t[i]), ground_truth)
}
points(1-speci, sensi, col = 'red', pch = 4)
lines(1-speci, sensi, col = 'red', lty = 3)