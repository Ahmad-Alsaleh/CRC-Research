"
This code uses clustering analysis to visualize if the identified features (OTUs)
can accurately cluster the patients into two distinct clusters: CRC and Normal.

"

# libraries ----
rm(list = ls())
library(stats)
library(ggplot2)
library(ggbiplot)
library(dplyr)
library(readxl)
source("Codes/Part 4 - Clustering Analysis/utility.R")

# PCA k-means clustering using all features ----

# reading and scaling the dataset
data = readxl::read_excel("Output Datasets/Analysis Dataset.xlsx") %>%
	(textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)

data[, -ncol(data)] = data[, -ncol(data)] %>% scale

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using All Features")
ggsave("Graphs/Clustering Using All Features.png", width = 2931, height = 1782, units = "px")

# PCA k-means clustering using important features ----

# extracting important features
features = read_excel("Output Datasets/Final Selected OTUs.xlsx")
features = features$OTU
data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using Significant Features")
ggsave("Graphs/Clustering Using Significant Features.png", width = 2931, height = 1782, units = "px")

