"
This code computes feature importance scores using WAB and BAM
framework. In this framework, different feature selection techniques,
like Information Gain, Symmetric Uncertainty, MRMR, Chi-Squared,
Random Forest, are aggregated to give a robust
list of the most important features.

After bootsrapping, the arithmetic mean is used to aggregate the
scores of each technique.

Importance scores are exported to: `Bootstrapped Importance Scores.xlsx`
and features names sorted based on importance scores are exported
to `Top Features.xlsx`.

This code is part of the CRC research.

Special thanks to Reem Salman.

Wednesday, 3 May 2023
Ahmad Alsaleh

"

# libraries ----
rm(list = ls())
source("Codes/Part 2 - Feature Selection/utility.R")
library(readxl)
library(splitstackshape)
library(FSelector)
library(praznik)
library(tibble)
library(xlsx)
library(caret)
library(dplyr)
library(ggplot2)
library(plotROC)
set.seed(42)

# loading data set ----
data = read_excel("Output Datasets/Analysis Dataset.xlsx") %>% 
	as.data.frame %>% (textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)
colnames(data)[ncol(data)] = "y"
str(data[, (ncol(data) - 2):ncol(data)])

# requirements for the wam() function ---

# Feature Selection Part 1: Filter methods ----
aggregated.scores = wam(data)


# TODO (Dr. Ayman) check if this is correct:
"
In WAM, scores of Random Forest is not included in WAB since it is not a filter
method but it is included in the BAM. (I have a feeling this is wrong).

"

# Feature Selection Part 2: Embedded and wrapper methods (machine learning) ----

# Random Forest ---
cv.folds = createFolds(data$y, k = 5, returnTrain = T)
training.control = trainControl(
	method = "cv",
	number = 5,
	search = "grid",
	classProbs = T,
	savePredictions = "final",
	index = cv.folds,
	summaryFunction = twoClassSummary
)

# fitting a random forest model with custom tuning parameters
rf.model = train(
	y ~ .,
	data = data,
	method = "rf",
	metric = "ROC",
	tuneGrid = expand.grid(.mtry = c(5, 25, 50, 100, 200, 300)),
	trControl = training.control,
	importance = T,
	nodesize = 1,
	ntree = 250,
	allowParallel = T
)
rf.model

# plotting ROC curve
roc.plot = ggplot(rf.model$pred[rf.model$pred$mtry == rf.model$finalModel$mtry,],
									aes(m = Normal, d = factor(obs, levels = c("Normal", "CRC")))) +
	geom_roc(n.cuts = 0) + coord_equal() + style_roc()
roc.plot = roc.plot + annotate("text", x = 0.75, y = 0.25, label =
															 	paste("AUC =", round((calc_auc(roc.plot))$AUC, 4)))
roc.plot

# plotting variables importance
var.imp.rf = varImp(rf.model)
plot(var.imp.rf, top = 20)

# extracting top variables from random forest model
var.imp.rf = var.imp.rf$importance[rownames(aggregated.scores),][, 1]
aggregated.scores = data.frame(aggregated.scores, `Random Forest` = var.imp.rf)

# BAM (aggregating using arithmetic mean) ----

# normalizing scores
aggregated.scores = as.data.frame(apply(aggregated.scores, 2, normalize))

# TODO (Dr. Ayman): is this how BAM should be computed? I am basically computing the 
# mean of all WAM scores. But note that I normalized before computing the mean.

# computing BAM
aggregated.scores$BAM = apply(aggregated.scores, 1, mean)

# sorting by score for each method ----
top.features = matrix(nrow = nrow(aggregated.scores),
											ncol = ncol(aggregated.scores))
for (col_i in 1:ncol(aggregated.scores))
	top.features[, col_i] = rownames(aggregated.scores[order(aggregated.scores[, col_i],
																													 decreasing = T),])

top.features = as.data.frame(top.features)
colnames(top.features) = c("IG", "SU", "MR", "CS", "RF", "BAM")

# exporting results ----
write.xlsx(aggregated.scores, "Output Datasets/Bootstrapped Importance Scores.xlsx")
write.xlsx(top.features, "Output Datasets/Top Features.xlsx", row.names = F)
