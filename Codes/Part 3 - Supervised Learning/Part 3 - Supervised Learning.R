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

Requirements:
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
data <- read_excel("Output Datasets/Analysis Dataset.xlsx") %>%
  (textshape::column_to_rownames)
colnames(data)[ncol(data)] <- "y"
data$y <- as.factor(data$y)

# reading the features file
imp.features.names <- read_excel("Output Datasets/Top Features.xlsx") %>% as.data.frame()

# fitting models and plotting AUC graphs ----

# SVM ---
SVM.predict.func <- function(model.fit, data.test) {
  attr(predict(model.fit, data.test, probability = T), "probabilities")[, levels(data.test$y)[2]]
}
SVM.AUC <- compute.AUCs(data, imp.features.names, SVM.predict.func, e1071::svm, probability = T)
SVM.plot <- generatePlot(SVM.AUC, "SVM")

# KNN (k = sqrt(n)) ---
KNN.predict.func <- function(model.fit, data.test) {
  predict(model.fit, data.test[, -ncol(data.test)], type = "prob")[, levels(data.test$y)[2]]
}
KNN.AUC <- compute.AUCs(data, imp.features.names, KNN.predict.func, e1071::gknn, k = floor(sqrt(nrow(data))))
KNN.plot <- generatePlot(KNN.AUC, expression(bold(paste("KNN (k = ", sqrt(n), ")"))))

# Naive Bayes ---
NB.predict.func <- function(model.fit, data.test) {
  predict(model.fit, data.test[, -ncol(data.test)], type = "raw")[, levels(data.test$y)[2]]
}
NB.AUC <- compute.AUCs(data, imp.features.names, NB.predict.func, e1071::naiveBayes)
NB.plot <- generatePlot(NB.AUC, "Naive Bayes")

# Random Forest ---
RF.predict.func <- function(model.fit, data.test) {
  predict(model.fit, data.test[, -ncol(data.test)], type = "prob")[, levels(data.test$y)[2]]
}
RF.AUC <- compute.AUCs(data, imp.features.names, RF.predict.func, randomForest::randomForest)
RF.plot <- generatePlot(RF.AUC, expression(bold(paste("Random Forest (m = ", sqrt(p), ")"))))

# combining all graphs in one grid ---
plot <- gridExtra::grid.arrange(SVM.plot, KNN.plot, NB.plot, RF.plot, nrow = 2, ncol = 2)
ggsave("Graphs/AUC Graphs.png", plot, width = 2931, height = 1782, units = "px", dpi = 210)

# extracting final subset of significant features ----
final.subset.name <- "BAM"

imp.features.names <- imp.features.names[, final.subset.name]

scores <- read_excel("Output Datasets/Bootstrapped Importance Scores.xlsx") %>%
  (textshape::column_to_rownames)
scores <- scores[imp.features.names, final.subset.name]

taxonomy <- read_excel("Output Datasets/taxonomy.xlsx") %>%
  (textshape::column_to_rownames)
taxonomy <- taxonomy[imp.features.names, c("Family", "Genus", "Species")]

final.subset <- data.frame(
  OTU = imp.features.names, `Imp Score` = scores,
  taxonomy, row.names = NULL
)

# using scree plot to choose cutoff for top features
scree.plot <- ggplot(
  data.frame(index = 1:nrow(final.subset), score = final.subset$Imp.Score),
  aes(x = index, y = score)
) +
  geom_line() +
  labs(x = "OTU Index", y = "Imp Score")
scree.plot

# drawing vertical line
cutoff <- 40

scree.plot + geom_vline(xintercept = cutoff, linetype = "dashed", color = "#D62828") +
  geom_text(
    x = cutoff, y = 0.01, hjust = -0.1, label = paste0("OTU Index = ", cutoff),
    color = "#D62828"
  )

ggsave("Graphs/BAM Scores Scree Plot.png", width = 2931, height = 891, units = "px")

# exporting final subset of significant features
final.subset <- final.subset[1:cutoff, ]
xlsx::write.xlsx(final.subset, "Output Datasets/Final Selected OTUs.xlsx", row.names = F)

# heatmap of relative abundance of the selected variables ----
data <- read_excel("Output Datasets/dataset.xlsx") %>%
  (textshape::column_to_rownames)

OTU <- t(data[, -(1:3)]) %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = T)
TAX <- read_excel("Output Datasets/taxonomy.xlsx") %>%
  (textshape::column_to_rownames) %>%
  as.matrix() %>%
  tax_table()

OTU <- OTU[final.subset$OTU, ]
TAX <- TAX[final.subset$OTU, ]

phy.seq <- phyloseq(OTU, TAX)
sample_data(phy.seq) <- data[, 1:3]

phy.seq <- prune_taxa(taxa_sums(phy.seq) > 0, phy.seq)
phy.seq <- transform_sample_counts(phy.seq, function(x) x / sum(x) * 100)

gpac <- prune_taxa(final.subset$OTU, phy.seq)
gplot <- plot_heatmap(gpac,
  sample.label = "Diagnosis", taxa.label = "Genus",
  sample.order = "Diagnosis", low = "#66CCFF", high = "#000033",
  na.value = "white"
)

gplot$labels$fill <- "Relative Abundance (%)"
gplot + geom_vline(xintercept = 47.5, color = "black", linewidth = 1.5, linetype = "dashed")

ggsave("Graphs/Relative Abundance of Significant Features.png", width = 2931, height = 1782, units = "px")
