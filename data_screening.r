##Breast Cancer Wisconsin Dataset
# Dua, D. and Graff, C. (2019). UCI Machine Learning Repository
# [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California,
# School of Information and Computer Science.

##load and clean up data
data <- read.csv("Datasets/breast-cancer-wisconsin.data", header = FALSE) #load dataset
colnames(data) <- c("ID", "Clump Thickness", "Uniformity of Cell Size",
                    "Uniformity of Cell Shape", "Marginal Adhesion",
                    "Single Epithelial Cell Size", "Bare Nuclei",
                    "Bland Chromatin", "Normal Nucleoli", "Mitoses",
                    "Diagnosis") #add descriptors
data$Diagnosis <- data$Diagnosis / 2 - 1 #normalize diagnosis label (0 = benign, 1 = malignant)
data <- data[!(data$`Bare Nuclei` == "?"),] #drop patients with missing observations
data$`Bare Nuclei` <- as.numeric(data$`Bare Nuclei`) #make all entries numeric

N <- length(data$Diagnosis)

##simulate 3 imperfect annotators
alpha <- c(0.99, 0.7, 0.6) #sensitivity
beta <- c(0.95, 0.7, 0.65) #specificity

Diagnosis1 <- numeric(length = N)
Diagnosis2 <- numeric(length = N)
Diagnosis3 <- numeric(length = N)

for(i in 1:N){
  if(data$Diagnosis[i] == 1){
    Diagnosis1[i] <- rbinom(1,1,alpha[1])
    Diagnosis2[i] <- rbinom(1,1,alpha[2])
    Diagnosis3[i] <- rbinom(1,1,alpha[3])
  } else {
    Diagnosis1[i] <- 1-rbinom(1,1,beta[1])
    Diagnosis2[i] <- 1-rbinom(1,1,beta[2])
    Diagnosis3[i] <- 1-rbinom(1,1,beta[3])
  }
}

#add to data
data$Diagnosis1 <- Diagnosis1
data$Diagnosis2 <- Diagnosis2
data$Diagnosis3 <- Diagnosis3

#create CSV file with simulated annotations
write.csv(data, "Datasets/BreastCancerWisconsinAnnotated06.csv")