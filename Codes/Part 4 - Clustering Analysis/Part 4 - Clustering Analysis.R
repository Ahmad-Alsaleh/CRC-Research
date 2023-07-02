"
This code uses clustering analysis to visualize if the identified features (OTUs)
can accurately cluster the patients into two distinct clusters: CRC and Normal.

"

# libraries ----
library(stats)
library(ggplot2)
library(ggbiplot)
library(dplyr)
library(readxl)
source("Codes/Part 4 - Clustering Analysis/utility.R")

# PCA k-means clustering using all features ----

# reading and scaling the dataset
data = readxl::read_excel(
	"Output Datasets/Analysis Dataset.xlsx") %>% 
	as.data.frame %>% (textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)

data[, -ncol(data)] = data[, -ncol(data)] %>% scale

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using All Features")

# PCA k-means clustering using important features ----

# extracting important features
features = read_excel("Output Datasets/Final Selected OTUs.xlsx")

# using scree plot to choose cutoff for top features
scree.plot = ggplot(data.frame(index = 1:nrow(features), score = features$score),
										aes(x = index, y = score)) + geom_line() + labs(x = "OTU Index", y = "Score")
scree.plot

cutoff = 35
scree.plot + geom_vline(xintercept = cutoff, linetype = "dashed", color = "red") +
	geom_text(x = cutoff, y = 0.1, hjust = -0.1, label = paste("Index =", cutoff), color = "red")

features = features$OTU[1:cutoff]
data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using Important Features")
