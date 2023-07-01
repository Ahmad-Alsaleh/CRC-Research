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
library(readxl)
library(splitstackshape)
library(FSelector)
library(praznik)
library(tibble)
library(xlsx)
library(caret)
library(dplyr)
set.seed(42)

# loading data set ----
data = read_excel("/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Analysis Dataset.xlsx") %>%
	as.data.frame
data = textshape::column_to_rownames(data)
data$Diagnosis = as.factor(data$Diagnosis)
colnames(data)[ncol(data)] = "y"
str(data[, (ncol(data) - 5):ncol(data)])

# helper functions ----
normalize = function(x) {
	return((x - min(x)) / (max(x) - min(x)))
}

# function that performs WAM on a single bootstrap ----

# feature selection techniques: Information Gain,
# Symmetric Uncertainty, MRMR, Chi-Squared

# Preconditions:
# The target variable is named "y"

wam = function(data) {
	iteration_i <<- iteration_i + 1
	
	# using different seed for each bootstrap
	set.seed(iteration_i)
	
	# getting stratified sample
	d = as.data.frame(stratified(data, "y", size = 0.999, replace = T))
	
	# Information Gain (IG)
	scores = information.gain(y ~ ., d)
	colnames(scores)[ncol(scores)] = "Information.Gain"
	
	# Symmetric Uncertainty (SU)
	scores = cbind(scores, symmetrical.uncertainty(y ~ ., d))
	colnames(scores)[ncol(scores)] = "Symmetric.Uncertainty"
	
	# Minimum Redundancy Maximum Relevance (MRMR)
	d.mrmr = MRMR(d[,!(colnames(d) %in% c("y"))], d$y, ncol(d) - 1)$score
	scores = merge(scores, d.mrmr, by = "row.names", all = T)
	scores = scores %>%
		remove_rownames() %>%
		column_to_rownames(var = "Row.names")
	colnames(scores)[ncol(scores)] = "MRMR"
	
	# Chi-Squared (CS)
	d.chi = chi.squared(y ~ ., d)
	scores = merge(scores, d.chi, by = "row.names", all = T)
	scores = scores %>%
		remove_rownames() %>%
		column_to_rownames(var = "Row.names")
	colnames(scores)[ncol(scores)] = "Chi.Squared"
	
	# updating the progress bar
	done_percent = round(iteration_i / num_bootstraps * 100, 1)
	text = paste("Iteration: ",
							 iteration_i,
							 "/",
							 num_bootstraps,
							 " (",
							 done_percent,
							 "%)",
							 sep = "")
	print(text)
	tcltk::setTkProgressBar(progress.bar, iteration_i, label = text)
	
	# return importance scores of a single bootstrap
	return(scores)
}

# requirements for the wam() function ---
num_bootstraps = 500
iteration_i = 0
progress.bar = tcltk::tkProgressBar(
	title = "Progress Bar",
	min = 0,
	max = num_bootstraps,
	width = 300
)

# Feature Selection Part 1: Filter methods ----
start.time = Sys.time()
bootstraps = replicate(num_bootstraps, wam(data), simplify = F)
bootstraps = abind::abind(bootstraps, along = 3)
excution.time = Sys.time() - start.time
excution.time

# aggregating results (WAM)
aggregated.scores = apply(bootstraps, c(1, 2), mean)

# TODO: check if this is correct:
"
In WAM, scores of Random Forest is not included in WAB since it is not a filter
method but it is included in the BAM. Is this correct? See the picture I sent.
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
library(ggplot2)
library(plotROC)
roc.plot = ggplot(rf.model$pred[rf.model$pred$mtry == rf.model$finalModel$mtry,],
									aes(m = Normal, d = factor(obs, levels = c("Normal", "CRC")))) +
	geom_roc(n.cuts = 0) + coord_equal() + style_roc()
roc.plot = roc.plot + annotate("text",
															 x = 0.75,
															 y = 0.25,
															 label = paste("AUC =", round((calc_auc(
															 	roc.plot
															 ))$AUC, 4)))
roc.plot

# plotting variables importance
var.imp.rf = varImp(rf.model)
plot(var.imp.rf, top = 20)

# extracting top variables from random forest model
var.imp.rf = var.imp.rf$importance[rownames(aggregated.scores),][, 1]
aggregated.scores = data.frame(aggregated.scores, `Random Forest` = var.imp.rf)

# BAM (aggregating using arithmetic mean) ----
aggregated.scores = as.data.frame(apply(aggregated.scores, 2, normalize))
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
write.csv(
	bootstraps,
	"/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Bootstraps.csv"
)
write.xlsx(
	aggregated.scores,
	"/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Bootstrapped Importance Scores.xlsx"
)
write.xlsx(
	top.features,
	"/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Top Features.xlsx",
	row.names = F
)
