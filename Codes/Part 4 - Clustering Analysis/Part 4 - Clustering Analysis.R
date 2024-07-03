"
This code uses clustering analysis to visualize if the identified features (OTUs)
can accurately cluster the patients into two distinct clusters: CRC and Normal.

"

# libraries ----
rm(list = ls())
library(stats)
library(ggplot2)
library(dplyr)
library(readxl)
source("Codes/Part 4 - Clustering Analysis/utility.R")

# PCA k-means clustering using all features ----

# reading and scaling the dataset
data <- readxl::read_excel("Output Datasets/Analysis Dataset.xlsx") %>%
  (textshape::column_to_rownames)
data$Diagnosis <- as.factor(data$Diagnosis)

data[, -ncol(data)] <- data[, -ncol(data)] %>% scale()

# clustering the observations
clusters.all.features <- kmeans.clusters(data, "PCA K-Means Clustering Using All Features")

# PCA k-means clustering using important features ----

# extracting important features
features <- read_excel("Output Datasets/Final Selected OTUs.xlsx")
features <- features$OTU
data <- data[, c(features, "Diagnosis")]

# clustering the observations
clusters.significant.features <- kmeans.clusters(data, "PCA K-Means Clustering Using Significant Features")

plot <- gridExtra::grid.arrange(
	clusters.all.features,
	clusters.significant.features,
  nrow = 2,
  ncol = 1
)
ggsave("Graphs/PCA K-Means Clustering.png", plot, width = 2931, height = 1782, units = "px")

