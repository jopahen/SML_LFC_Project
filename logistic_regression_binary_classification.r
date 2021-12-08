set.seed(10)
#load data, we work with R = 3 simulated erroneous annotators, 80/20 train/test split
#data = training data, data_test = testing data
data <- read.csv("Datasets/BreastCancerWisconsinAnnotated07normalized.csv")[,-1]
inds <- sample(1:length(data$Diagnosis), size = floor(0.8 * length(data$Diagnosis)),
               replace = FALSE)
data_test <- data[-inds,]
data <- data[inds,]
R <- 3

options(warn = -1)
#prepare data for ground truth logistic regression, fit logistic regression and display predictions
data_ground_truth <- data[,-c(1, 12, 13, 14)]
model_ground_truth <- glm(Diagnosis ~ . - 1, family = binomial, data = data_ground_truth)
as.numeric(model_ground_truth$fitted.values > 0.5)

#prepare data for majority vote logistic regression, fit logistic regression and display predictions
Ann <- as.matrix(t(cbind(data$Diagnosis1,data$Diagnosis2,data$Diagnosis3)))
Diagnosis_majority <- as.numeric(colSums(Ann) / R > 0.5)
Diagnosis_majority_probs <- colSums(Ann) / R
data_majority <- data[,-c(1, 11, 12, 13, 14)]
data_majority$DiagnosisMajority <- Diagnosis_majority
model_majority <- glm(DiagnosisMajority ~ . - 1, family = binomial, data = data_majority)
as.numeric(model_majority$fitted.values > 0.5)
w_init_log <- as.vector(model_majority$coefficients)
colnames(w_init_log) <- NULL

options(warn = 0)
