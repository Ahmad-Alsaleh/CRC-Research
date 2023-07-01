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
				 subtitle = sprintf("Accuracy: %.1f%%", find.clustering.accuracy(clustering.scores) * 100)) +
		geom_point(size = 5, alpha = 0.8) + ggtitle(title) +
		theme(plot.title = element_text(hjust = 0.5)) +
		geom_point(size = 1, color = "black", alpha = 0.2)
}

# function to run K-means clustering on PCA scores and then plotting
# the clusters in 2D using the first two clusters
kmeans.clusters = function(data, title) {
	# finding PCA scores
	pca.scores = prcomp(data[, -ncol(data)])$x %>% as.data.frame
	
	# k-means clustering
	set.seed(42)
	clusters.labels = kmeans(pca.scores, centers = 2, nstart = 20)$cluster
	pca.scores$Cluster = as.factor(clusters.labels)
	pca.scores$Diagnosis = data$Diagnosis
	
	# plotting clusters
	plot.clusters(pca.scores, title)
	
}
