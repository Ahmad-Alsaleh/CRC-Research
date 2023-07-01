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
	"/Users/ahmad/AUS/Research/CRC Research/Output Datasets/Analysis Dataset.xlsx") %>% 
	as.data.frame %>% (textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)

# TODO: try without this
data[, -ncol(data)] = data[, -ncol(data)] %>% scale

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using All Features")

# PCA k-means clustering using important features ----

# extracting important features
features = read_excel("/Users/ahmad/AUS/Research/CRC Research/Output Datasets/Final Selected OTUs.xlsx")

# using scree plot to choose cutoff for top features
features = features$OTU[1:30] # TODO: make this as 60 when exporting the file (use a scree plot with a vertical line)

data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using Important Features")
