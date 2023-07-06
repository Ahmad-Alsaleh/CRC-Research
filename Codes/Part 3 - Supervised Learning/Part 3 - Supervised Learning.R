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
"
-	Response must be the last column of the dataset
-	The columns of the features file must contain the different features
	names ranked from the most to the least important.

"

# libraries ----
set.seed(42)
rm(list = ls())
library(dplyr)
library(phyloseq)
library(readxl)
source("Codes/Part 3 - Supervised Learning/utility.R")

# reading files ----

# reading the dataset file
data = read_excel("Output Datasets/Analysis Dataset.xlsx") %>% as.data.frame %>%
	(textshape::column_to_rownames)
colnames(data)[ncol(data)] = "y"
data$y = as.factor(data$y)

# reading the features file
imp.features = read_excel("Output Datasets/Top Features.xlsx") %>% as.data.frame

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
taxonomy = read_excel("Output Datasets/taxonomy.xlsx") %>% 
	(textshape::column_to_rownames)

scores = read_excel("Output Datasets/Bootstrapped Importance Scores.xlsx") %>% 
	(textshape::column_to_rownames) %>% as.data.frame

scores = scores[final.subset, final.subset.name]

final.featrues = data.frame(OTU = final.subset, `BAM Score` = scores,
														taxonomy[final.subset, c("Family", "Genus", "Species")])
xlsx::write.xlsx(final.featrues, "Output Datasets/Final Selected OTUs.xlsx", row.names = F)

# heatmap of relative abundance of the selected variables
data = read_excel("Output Datasets/dataset.xlsx") %>% (textshape::column_to_rownames)

OTU = t(data[, -(1:3)]) %>% as.matrix %>% otu_table(taxa_are_rows = T)
TAX = read_excel("Output Datasets/taxonomy.xlsx") %>% (textshape::column_to_rownames) %>%
	as.matrix %>% tax_table

phy.seq = phyloseq(OTU, TAX)
sample_data(phy.seq) = data[, 1:3]

phy.seq = prune_taxa(taxa_sums(phy.seq) > 0, phy.seq)
gpac = prune_taxa(final.featrues$OTU[1:50], phy.seq)

gplot = plot_heatmap(gpac, sample.label = "Diagnosis",taxa.label = "Genus",
										 sample.order = "Diagnosis", low="#66CCFF", high="#000033",
										 na.value="white")

gplot$labels$fill = "Relative Abundance (%)"
gplot + geom_vline(xintercept = 47.5, color = "black", linewidth = 1.5, linetype = "dashed") +
	ggtitle("Relative Abundance of Top 50 Features")
