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
scree.plot = ggplot(data.frame(index = 1:nrow(features), score = features$BAM.Score),
										aes(x = index, y = score)) + geom_line() + labs(x = "OTU Index", y = "BAM Score")
scree.plot

# horizontal line
h.cutoff = 0.22
scree.plot + geom_hline(yintercept = h.cutoff, linetype = "dashed", color = "red")

v.cutoff = which(features$BAM.Score < h.cutoff)[1]

# vertical line
scree.plot + geom_vline(xintercept = v.cutoff, linetype = "dashed", color = "red") +
	geom_hline(yintercept = h.cutoff, linetype = "dashed", color = "red") +
	geom_text(x = v.cutoff, y = 0.1, hjust = -0.1, label = paste("OTU Index =", v.cutoff), color = "red")

features = features$OTU[1:v.cutoff]
data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data, "PCA K-Means Clustering Using Important Features")
