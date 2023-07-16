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
library(xlsx)
library(dplyr)
set.seed(42)

# loading data set ----
data = read_excel("Output Datasets/Analysis Dataset.xlsx") %>%
	as.data.frame %>% (textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)
colnames(data)[ncol(data)] = "y"
str(data[, (ncol(data) - 2):ncol(data)])

# feature selection ----
aggregated.scores = aggregate.scores(data)

# sorting OTUs by importance scores of each method ----
top.features = matrix(nrow = nrow(aggregated.scores), ncol = ncol(aggregated.scores))
colnames(top.features) = colnames(aggregated.scores)

for (col_i in 1:ncol(aggregated.scores))
	top.features[, col_i] = rownames(aggregated.scores)[order(-aggregated.scores[, col_i])]

# exporting results ----
write.xlsx(aggregated.scores, "Output Datasets/Bootstrapped Importance Scores.xlsx")
write.xlsx(top.features, "Output Datasets/Top Features.xlsx", row.names = F)

