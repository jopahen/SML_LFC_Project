# SML_LFC_Project
Code for implementations of the SML-Project on "Learning from Crowds" on the paper by Raykar et al.

run data_screening.r first to simulate fresh annotations + (un-)normalized data files (.csv)

running raykar_binary_classification.r fits the model using EM-algorithm and displays estimates/predictions/misclassifications at the end

running logistic_regression_binary_classification.r fits simple logistic regression models on majority + ground truth and displays predictions/misclassifications at the end

running performance.r will run a pipeline that fits raykar and simple majority vote logistic model and evaluates performance using ROC-plots, gives out confusion matrices, performance on test data during training

running run.r runs everything from start to finish
