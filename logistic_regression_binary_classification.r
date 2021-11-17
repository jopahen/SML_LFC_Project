#load data, we work with R = 3 simulated erroneous annotators
data <- read.csv("Datasets/BreastCancerWisconsinAnnotated.csv")[,-1]
R <- 3

#prepare data for ground truth logistic regression, fit logistic regression and display predictions
data_ground_truth <- data[,-c(1, 12, 13, 14)]
model_ground_truth <- glm(Diagnosis ~ ., family = binomial, data = data_ground_truth)
as.numeric(model_ground_truth$fitted.values > 0.5)

#prepare data for majority vote logistic regression, fit logistic regression and display predictions
Ann <- as.matrix(t(cbind(data$Diagnosis1,data$Diagnosis2,data$Diagnosis3)))
Diagnosis_majority <- colSums(Ann) / R
data_majority <- data[,-c(1, 11, 12, 13, 14)]
data_majority$DiagnosisMajority <- Diagnosis_majority
model_majority <- glm(DiagnosisMajority ~ ., family = binomial, data = data_majority)
as.numeric(model_majority$fitted.values > 0.5)

#eyeball misclassification error
sum(data$Diagnosis != as.numeric(model_ground_truth$fitted.values > 0.5))
sum(data$Diagnosis != as.numeric(model_majority$fitted.values > 0.5))
