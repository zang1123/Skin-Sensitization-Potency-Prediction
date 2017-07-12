# Open R and set the working environment to "C:/SkinSensitizationClassification".

setwd("C:/SkinSensitizationClassification")

# Install and load package e1071 which contains the svm function.

install.packages("e1071")
library(e1071)

# I. LLNA Classification

# Read in the LLNA data with 120 chemicals and 16 columns.

LLNA <- read.table("LLNA-SkinData.txt", header=T, sep="\t", as.is=T )

# LLNA Sensitizers vs Non-sensitizers Model. 

# Use eight variables (hCLAT, OECD and six properties) to build models.

# There are 94 chemicals in training set (rows 1-94).  

LLNASensNonsensTraining <- LLNA[1:94, c("hCLAT", "OECD", "MW", "LogP", "LogS", "LogVP", "MP", "BP")]

# There are 26 chemicals in test set (rows 95-120). 

LLNASensNonsensTest <- LLNA[95:120, c("hCLAT", "OECD", "MW", "LogP", "LogS", "LogVP", "MP", "BP")]

# There are 68 positives (sensitizers) and 26 negatives (non-sensitizers) in training set

TrainingClassSensNonsens <- factor(c(rep("POS", 68), rep("NEG", 26)))

# There are 19 positives (sensitizers) and 7 (non-sensitizers) negatives in test set

TestClassSensNonsens <- factor(c(rep("POS", 19), rep("NEG", 7)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value from a series of trials from 1 to 10.

wtsSensNonsens <- 7/ table(TrainingClassSensNonsens)

# Build support vector machine model using svm function. Two parameters cost and gamma are set to 4600 and 0.022, respectively via an optimization procedure.

SVMmodelSensNonsens <- svm(LLNASensNonsensTraining, TrainingClassSensNonsens, cost = 4600, gamma = 0.022, class.weights = wtsSensNonsens)

# Predict the training set and produce a confusion matrix. 

PredTrainSensNonsens <-predict(SVMmodelSensNonsens, LLNASensNonsensTraining)

table(PredTrainSensNonsens, TrainingClassSensNonsens)

# Predict the test set and produce a confusion matrix. 

PredTestSensNonsens <-predict(SVMmodelSensNonsens, LLNASensNonsensTest)

table(PredTestSensNonsens, TestClassSensNonsens)

# Leave-one-out cross validation. Use 94 training chemicals (set the parameter cross=94). 

SVMmodelSensNonsensLOOCV <- svm(LLNASensNonsensTraining, TrainingClassSensNonsens, cross = 94, cost = 4600, gamma = 0.022, class.weights = wtsSensNonsens)

summary(SVMmodelSensNonsensLOOCV)


# LLNA Strong potency vs weak potency Model.

# Use nine variables (three assays and six properties) to build models.

# There are 68 chemicals in training set (rows 1-68).  

LLNAStrongWeakTraining <- LLNA[1:68, c("MW", "LogP", "LogS", "LogVP", "MP", "BP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 19 chemicals in test set (rows 95-113). 

LLNAStrongWeakTest <- LLNA[95:113, c("MW", "LogP", "LogS", "LogVP", "MP", "BP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 28 1A (strong sensitizers) and 40 1B (weak sensitizers) in training set.

TrainingClassStrongWeak <- factor(c(rep("A", 28), rep("B", 40)))

# There are 7 1A (strong sensitizers) and 12 1B (weak sensitizers) negatives in test set

TestClassStrongWeak <- factor(c(rep("A", 7), rep("B", 12)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value #from a series of trials from 1 to 10.

wtsStrongWeak <- 7/ table(TrainingClassStrongWeak)

# Build support vector machine model using svm function. Two parameters cost and gamma are #set to 800 and 0.002, respectively via an optimization procedure.

SVMmodelStrongWeak <- svm(LLNAStrongWeakTraining, TrainingClassStrongWeak, cost = 800, gamma = 0.002, class.weights = wtsStrongWeak)

# Predict the training set and produce a confusion matrix. 

PredTrainStrongWeak <-predict(SVMmodelStrongWeak, LLNAStrongWeakTraining)

table(PredTrainStrongWeak, TrainingClassStrongWeak)

# Predict the test set and produce a confusion matrix. 

PredTestStrongWeak <-predict(SVMmodelStrongWeak, LLNAStrongWeakTest)

table(PredTestStrongWeak, TestClassStrongWeak)

# Leave-one-out cross validation. Use 68 training chemicals (set the parameter cross=68). 

SVMmodelStrongWeakLOOCV <- svm(LLNAStrongWeakTraining, TrainingClassStrongWeak, cross = 68, cost = 800, gamma = 0.002, class.weights = wtsStrongWeak)

summary(SVMmodelStrongWeakLOOCV)


# LLNA Three-Category Classification Model.  

# Use nine variables (three assays and six properties) to build models.

# There are 94 chemicals in training set (rows 1-94). 

LLNAThreeClassTraining <- LLNA[1:94, c("MW", "LogP", "LogS", "LogVP", "MP", "BP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 26 chemicals in test set (rows 95-120). 

LLNAThreeClassTest <- LLNA[95:120, c("MW", "LogP", "LogS", "LogVP", "MP", "BP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 26 negatives and 68 positives (1A=28 and 1B=40) for training set

TrainingClassThreeClass <- factor(c(rep("A", 28), rep("B", 40), rep("NEG", 26)))

# There are 7 negatives and 19 positives (A=7 and B=12) for test set

TestClassThreeClass <- factor(c(rep("A", 7), rep("B", 12), rep("NEG", 7)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value #from a series of trials from 1 to 10.

wtsThreeClass <- 7/ table(TrainingClassThreeClass)

# Build support vector machine model using svm function. Two parameters cost and gamma are #set to 1500 and 0.01, respectively via an optimization procedure.

SVMmodelThreeClass <- svm(LLNAThreeClassTraining, TrainingClassThreeClass, cost = 1500, gamma = 0.01, class.weights = wtsThreeClass)

# Predict the training set and produce a confusion matrix. 

PredTrainThreeClass <-predict(SVMmodelThreeClass, LLNAThreeClassTraining)

table(PredTrainThreeClass, TrainingClassThreeClass)

# Predict the test set and produce a confusion matrix. 

PredTestThreeClass <-predict(SVMmodelThreeClass, LLNAThreeClassTest)

table(PredTestThreeClass, TestClassThreeClass)

# Leave-one-out cross validation. Use 94 training chemicals (set the parameter cross=94). 

SVMmodelThreeClassLOOCV <- svm(LLNAThreeClassTraining, TrainingClassThreeClass, cross = 94, cost = 1500, gamma = 0.01, class.weights = wtsThreeClass)

summary(SVMmodelThreeClassLOOCV)


# II. Human Classification

# Human Data for Sensitizers vs Non-sensitizers. Read in the human data with 96 chemicals and 17 columns.

HumanSensNonsens <- read.table("Human-Sensitizers-Nonsensitizers.txt", header=T, sep="\t", as.is=T )


# Human Data for strong potency vs weak potency and three category modeling. Read in the human data with 87 chemicals and 14 columns.

HumanThreeCategories <- read.table("Human-ThreeCategories.txt", header=T, sep="\t", as.is=T )


# Human Sensitizers vs Non-sensitizers Model

# Use five variables (hCLAT, DPRA, KeratinoSens, OECD and LogP) to build models.

# There are 72 chemicals in training set (rows 1-72).  

HumanSensNonsensTraining <- HumanSensNonsens[1:72, c("LogP", "avg.Lys.Cys",  "hCLAT", "Keratino", "OECD")]

# There are 24 chemicals in test set (rows 73-96). 

HumanSensNonsensTest <- HumanSensNonsens[73:96, c("LogP", "avg.Lys.Cys",  "hCLAT", "Keratino", "OECD")]

# There are 21 negatives (non-sensitizers) and 51 positives (sensitizers) in training set.

TrainingClassSensNonsens <- factor(c(rep("NEG", 21), rep("POS", 51)))

# There are 9 negatives (non-sensitizers) and 15 positives (sensitizers) in test set.

TestClassSensNonsens <- factor(c(rep("NEG", 9), rep("POS", 15)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value #from a series of trials from 1 to 10.

wtsSensNonsens <- 7/ table(TrainingClassSensNonsens)

# Build support vector machine model using svm function. Two parameters cost and gamma are #set to 1000 and 0.05, respectively via an optimization procedure.

SVMmodelSensNonsens <- svm(HumanSensNonsensTraining, TrainingClassSensNonsens, cost = 1000, gamma = 0.05, class.weights = wtsSensNonsens)

# Predict the training set and produce a confusion matrix. 

PredTrainSensNonsens <-predict(SVMmodelSensNonsens, HumanSensNonsensTraining)

table(PredTrainSensNonsens, TrainingClassSensNonsens)

# Predict the test set and produce a confusion matrix. 

PredTestSensNonsens <-predict(SVMmodelSensNonsens, HumanSensNonsensTest)

table(PredTestSensNonsens, TestClassSensNonsens)

# Leave-one-out cross validation. Use 72 training chemicals (set the parameter cross=72). 

SVMmodelSensNonsensLOOCV <- svm(HumanSensNonsensTraining, TrainingClassSensNonsens, cross = 72, cost = 1000, gamma = 0.05, class.weights = wtsSensNonsens)

summary(SVMmodelSensNonsensLOOCV)


# Human Strong potency vs weak potency Model

# Use four variables (three assays and LogP) to build models.

# There are 41 chemicals in training set (rows 1-41).  

HumanStrongWeakTraining <- HumanThreeCategories[1:41, c("LogP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 16 chemicals in test set (rows 64-79). 

HumanStrongWeakTest <- HumanThreeCategories[64:79, c("LogP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 19 1A (strong sensitizers) and 22 1B (weak sensitizers) in training set.

TrainingClassStrongWeak <- factor(c(rep("A", 19), rep("B", 22)))

# There are 7 1A (strong sensitizers) and 9 1B (weak sensitizers) negatives in test set.

TestClassStrongWeak <- factor(c(rep("A", 7), rep("B", 9)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value #from a series of trials from 1 to 10.

wtsStrongWeak <- 7/ table(TrainingClassStrongWeak)

# Build support vector machine model using svm function. Two parameters cost and gamma are #set to 1000 and 0.016, respectively via an optimization procedure.

SVMmodelStrongWeak <- svm(HumanStrongWeakTraining, TrainingClassStrongWeak, cost = 1000, gamma = 0.016, class.weights = wtsStrongWeak)

# Predict the training set and produce a confusion matrix. 

PredTrainStrongWeak <-predict(SVMmodelStrongWeak, HumanStrongWeakTraining)

table(PredTrainStrongWeak, TrainingClassStrongWeak)

# Predict the test set and produce a confusion matrix. 

PredTestStrongWeak <-predict(SVMmodelStrongWeak, HumanStrongWeakTest)

table(PredTestStrongWeak, TestClassStrongWeak)

# Leave-one-out cross validation. Use 41 training chemicals (set the parameter cross=41). 

SVMmodelStrongWeakLOOCV <- svm(HumanStrongWeakTraining, TrainingClassStrongWeak, cross = 41, cost = 1000, gamma = 0.016, class.weights = wtsStrongWeak)

summary(SVMmodelStrongWeakLOOCV)


# Human Three-Category Classification Model.  

# Use four variables (three assays and LogP) to build models.

# There are 63 chemicals in training set (rows 1-63). 

HumanThreeClassTraining <- HumanThreeCategories[1:63, c("LogP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 24 chemicals in test set (rows 64-87). 

HumanThreeClassTest <- HumanThreeCategories[64:87, c("LogP", "hCLAT.CD86.CD54.MIN", "avgDepletionLysCys", "ARE.EC1.5")]

# There are 22 negatives and 41 positives (1A=19 and 1B=22) for training set.

TrainingClassThreeClass <- factor(c(rep("A", 19), rep("B", 22), rep("NEG", 22)))

# There are 8 negatives and 16 positives (A=7 and B=9) for test set.

TestClassThreeClass <- factor(c(rep("A", 7), rep("B", 9), rep("NEG", 8)))

# Set a weight for balancing the negative and positive classes, where "7" is an optimal value #from a series of trials from 1 to 10.

wtsThreeClass <- 7/ table(TrainingClassThreeClass)

# Build support vector machine model using svm function. Two parameters cost and gamma are #set to 80 and 0.011, respectively via an optimization procedure.

SVMmodelThreeClass <- svm(HumanThreeClassTraining, TrainingClassThreeClass, cost = 80, gamma = 0.011, class.weights = wtsThreeClass)

# Predict the training set and produce a confusion matrix. 

PredTrainThreeClass <-predict(SVMmodelThreeClass, HumanThreeClassTraining)

table(PredTrainThreeClass, TrainingClassThreeClass)

# Predict the test set and produce a confusion matrix. 

PredTestThreeClass <-predict(SVMmodelThreeClass, HumanThreeClassTest)

table(PredTestThreeClass, TestClassThreeClass)

# Leave-one-out cross validation. Use 63 training chemicals (set the parameter cross=63). 

SVMmodelThreeClassLOOCV <- svm(HumanThreeClassTraining, TrainingClassThreeClass, cross = 63, cost = 80, gamma = 0.011, class.weights = wtsThreeClass)

summary(SVMmodelThreeClassLOOCV)





