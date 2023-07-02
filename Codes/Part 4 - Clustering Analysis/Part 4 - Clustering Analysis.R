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
# TODO: in Part 3, add a column of the scores values (mostly make it col 2)
features = read_excel("Output Datasets/Final Selected OTUs.xlsx")

# using scree plot to choose cutoff for top features

# TODO: (use a scree plot with a vertical line).
cutoff = 30
features = features$OTU[1:cutoff]

data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using Important Features")
