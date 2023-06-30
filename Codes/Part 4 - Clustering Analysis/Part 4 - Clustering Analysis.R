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

# functions ----

# function to calculate the percent of observations in the same cluster
find.clustering.accuracy = function(clustering.scores) {
	t = table(clustering.scores$Cluster, clustering.scores$Diagnosis)
	(sum(diag(t)) / sum(t)) %>% 
		ifelse(. < 0.5, 1 - ., .)
}

# function to generate a plot with predefined custom aesthetics
plot.clusters = function(clustering.scores, title) {
	ggplot(clustering.scores, aes(x = PC1, y = PC2, color = Cluster, shape = Diagnosis)) +
		labs(x = "Principal Component 1", y = "Principal Component 2",
				 subtitle = sprintf("Accuracy: %.2f", find.clustering.accuracy(clustering.scores))) +
		geom_point(size = 3) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}

# function to run K-means clustering on PCA scores and then plotting
# the clusters in 2D using the first two clusters
kmeans.clusters = function(data) {
	# finding PCA scores
	pca.scores = prcomp(data[, -ncol(data)])$x %>% as.data.frame
	
	# k-means clustering
	set.seed(42)
	clusters.labels = kmeans(pca.scores, centers = 2, nstart = 20)$cluster
	pca.scores$Cluster = as.factor(clusters.labels)
	pca.scores$Diagnosis = data$Diagnosis
	
	# plotting clusters
	plot.clusters(pca.scores, "PCA K-Means Clustering Using Important Features")
	
}

# reading dataset ----
data = readxl::read_excel(
	"/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Analysis Dataset.xlsx") %>% 
	as.data.frame %>% (textshape::column_to_rownames)
data$Diagnosis = as.factor(data$Diagnosis)

data[, -ncol(data)] = data[, -ncol(data)] %>% scale

# PCA k-means clustering using all features ----

# clustering the observations
kmeans.clusters(data)


# PCA k-means clustering using important features ----

# extracting important features
features = read_excel("/Users/ahmad/AUS/Research/Cancer Project/Output Datasets/Final Selected OTUs.xlsx")
features = features$OTU[1:30] # TODO: make this as 60 when exporting the file (use a scree plot with a vertical line)

data = data[, c(features, "Diagnosis")]

# clustering the observations
kmeans.clusters(data)

# ----


























