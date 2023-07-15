# function that performs WAM on a single bootstrap ----

# feature selection techniques: Information Gain,
# Symmetric Uncertainty, MRMR, Chi-Squared

# Preconditions:
# - The target variable is named "y"
# - `iteration_i` must increase in each iteration as it is used for setting the seed
single.bootstrap_ = function(data, iteration_i, num_bootstraps) {
	# using different seed for each bootstrap
	set.seed(iteration_i)
	
	# getting stratified sample
	d = stratified(data, "y", size = 0.999, replace = T)
	
	# importance scores of the features
	scores = matrix(nrow = ncol(d) - 1, ncol = 5, dimnames =
										list(colnames(d)[1:(ncol(d)-1)], c("IG", "SU", "GR", "CS", "BAM")))
	
	# Information Gain (IG)
	scores[, "IG"] = information.gain(y ~ ., d)[, 1]
	
	# Symmetric Uncertainty (SU)
	scores[, "SU"] = symmetrical.uncertainty(y ~ ., d)[, 1]
	
	# Gain Ratio (GR)
	scores[, "GR"] = gain.ratio(y ~ ., d)[, 1]
	
	# Chi-Squared (CS)
	scores[, "CS"] = chi.squared(y ~ ., d)[, 1]
	
	# normalizing scores of each feature selection method
	scores[, 1:4] = apply(scores[, 1:4], 2, normalize)
	
	# BAM
	scores[, "BAM"] = apply(scores[, 1:4], 1, mean)
	
	# printing updated progress
	done_percent = round(iteration_i / num_bootstraps * 100, 1)
	text = paste0("Iteration: ", iteration_i, "/", num_bootstraps, " (",
								done_percent, "%)")
	print(text)
	
	# return importance scores of a single bootstrap
	return(scores)
}

# normalizes a vector using min-max normalization ----
normalize = function(x) {
	return((x - min(x)) / (max(x) - min(x)))
}

# performs WAM on `data` ----
# Preconditions:
# - The target variable should be the last column in `data` and should be named "y"
aggregate.scores = function(data, num_bootstraps = 500) {
	iteration_i = 0
	
	# bootstrapping `num_bootstraps` times
	bootstraps = replicate(num_bootstraps, simplify = F, expr = { # TODO (Ahmad Alsaleh): try to use multi-core processing for replicate()
		iteration_i <<- iteration_i + 1
		single.bootstrap_(data, iteration_i, num_bootstraps)
	}) %>% abind::abind(along = 3)
	
	# aggregating scores of all bootstraps using the arithmetic mean
	aggregated.bootstraps = apply(bootstraps, 1:2, mean)
	
	return(aggregated.bootstraps)
}
