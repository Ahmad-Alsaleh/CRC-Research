# Code description ----
"
This code runs differnt supervised machine learning models on the dataset
using differnt subsets of features to assess the accuracy of prediction
using AUC as a measure.

The ML models used in this code are: KNN (k = sqrt(n)), Naive Bayes,
Random Forest, and SVM.

This code is part of the CRC research paper.

Special thanks to Reem Salman.

Ahmad Alsaleh
Under the supervision of Dr. Ayman Alzaatreh.
May 27th 2022.

"

# Requirements ----
# TODO: remove this (or some of it)
"
-	Change the paths of the dataset and the features files (both files shluld be
	of type xlsx).
-	Please note that the response must be the last column of the dataset
-	The columns of the features excel file contains the different features
	names ranked from the most to the least important.
-	(Optional) modify the values of the number of CV folds and the number of
	points to show on the plot. Be aware that increasing any of these two values
	might slow the code dramatically and it is thus not recommended to change
	the default values unless needed.

"

# libraries ----
set.seed(42)
library(dplyr)
source("Codes/Part 3 - Supervised Learning/utility.R")

# reading files ----

# reading the dataset file
data = readxl::read_excel("Output Datasets/Analysis Dataset.xlsx") %>% as.data.frame %>%
	(textshape::column_to_rownames)
colnames(data)[ncol(data)] = "y"
data$y = as.factor(data$y)

# reading the features file
imp.features = readxl::read_excel("Output Datasets/Top Features.xlsx") %>% as.data.frame

# fitting models and plotting AUC graphs ----

# SVM ---
SVM.predict.func = function(model.fit, data.test)
	attr(predict(model.fit, data.test, probability = T), "probabilities")[, levels(data.test$y)[2]]
SVM.AUC = compute.AUCs(data, imp.features, SVM.predict.func, e1071::svm, probability = T)
SVM.plot = generatePlot(SVM.AUC, "SVM")

# KNN (k = sqrt(n)) ---
KNN.predict.func = function(model.fit, data.test)
	# TODO: check if [1] or [2] is needed bellow (print it in the function). [1] gives better graph
	predict(model.fit, data.test[, -ncol(data.test)], type = "prob")[, levels(data.test$y)[1]]
KNN.AUC = compute.AUCs(data, imp.features, KNN.predict.func, e1071::gknn, k = floor(sqrt(nrow(data))))
KNN.plot = generatePlot(KNN.AUC, expression(bold(paste("KNN (k = ", sqrt(n), ")"))))

# Naive Bayes ---
NB.predict.func = function(model.fit, data.test)
	predict(model.fit, data.test[, -ncol(data.test)], type = "raw")[, levels(data.test$y)[2]]
NB.AUC = compute.AUCs(data, imp.features, NB.predict.func, e1071::naiveBayes)
NB.plot = generatePlot(NB.AUC, "Naive Bayes")

# Random Forest ---
RF.predict.func = function(model.fit, data.test)
	predict(model.fit, data.test[, -ncol(data.test)], type = "prob")[, levels(data.test$y)[2]]
RF.AUC = compute.AUCs(data, imp.features, RF.predict.func, randomForest::randomForest)
RF.plot = generatePlot(RF.AUC, expression(bold(paste("Random Forest (m = ", sqrt(p), ")"))))

# combining all graphs in one grid ---
gridExtra::grid.arrange(SVM.plot, KNN.plot, NB.plot, RF.plot, nrow = 2, ncol = 2)

# extracting final features ----
final.subset.name = "BAM"

final.subset = imp.features[, final.subset.name]
taxonomy = readxl::read_excel("Output Datasets/taxonomy.xlsx") %>% 
	(textshape::column_to_rownames)

scores = readxl::read_excel("Output Datasets/Bootstrapped Importance Scores.xlsx") %>% 
	(textshape::column_to_rownames) %>% as.data.frame
scores = scores[final.subset, final.subset.name]

final.featrues = data.frame(OTU = final.subset, score = scores,
														taxonomy[final.subset, c("Family", "Genus", "Species")])
xlsx::write.xlsx(final.featrues, "Output Datasets/Final Selected OTUs.xlsx", row.names = F)

# heatmap of relative abundance of the selected variables
# TODO:
"
re-create `phy.seq` object here. this snippet of code requires the `phy.seq`
object from Part 1

"
gpac = prune_taxa(final.featrues$OTU[1:50], phy.seq)
gplot = plot_heatmap(gpac, sample.label = "Diagnosis",taxa.label = "Genus",
										 sample.order = "Diagnosis", low="#66CCFF", high="#000033",
										 na.value="white")

gplot$labels$fill = "Relative Abundance (%)"
gplot + geom_vline(xintercept = 47.5, color = "black", size = 1.5, linetype = "dashed") +
	ggtitle("Relative Abundance of Top 50 Features")

