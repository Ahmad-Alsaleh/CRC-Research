# function that performs WAM on a single bootstrap ----

# feature selection techniques: Information Gain,
# Symmetric Uncertainty, MRMR, Chi-Squared

# Preconditions:
# - The target variable is named "y"
# - `iteration_i` must increase in each iteration as it is used for setting the seed
single.wam_ = function(data, iteration_i, num_bootstraps) {
	# using different seed for each bootstrap
	set.seed(iteration_i)
	
	# getting stratified sample
	d = as.data.frame(stratified(data, "y", size = 0.999, replace = T))
	
	# Information Gain (IG)
	scores = information.gain(y ~ ., d)
	colnames(scores)[ncol(scores)] = "Information.Gain"
	
	# Symmetric Uncertainty (SU)
	scores = cbind(scores, symmetrical.uncertainty(y ~ ., d))
	colnames(scores)[ncol(scores)] = "Symmetric.Uncertainty"
	
	# Minimum Redundancy Maximum Relevance (MRMR)
	d.mrmr = MRMR(d[,!(colnames(d) %in% c("y"))], d$y, ncol(d) - 1)$score
	scores = merge(scores, d.mrmr, by = "row.names", all = T)
	scores = scores %>%
		remove_rownames() %>%
		column_to_rownames(var = "Row.names")
	colnames(scores)[ncol(scores)] = "MRMR"
	
	# Chi-Squared (CS)
	d.chi = chi.squared(y ~ ., d)
	scores = merge(scores, d.chi, by = "row.names", all = T)
	scores = scores %>%
		remove_rownames() %>%
		column_to_rownames(var = "Row.names")
	colnames(scores)[ncol(scores)] = "Chi.Squared"
	
	# printing updated progress
	done_percent = round(iteration_i / num_bootstraps * 100, 1)
	text = paste("Iteration: ", iteration_i, "/", num_bootstraps, " (",
							 done_percent, "%)", sep = "")
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
# - The target variable is named "y"
# TODO (Dr. Ayman):
# Please check if WAM is implemented correctly, especially if
# the aggreagation and normalization are done in the right place.
# I am using `apply(bootstraps, 1:2, mean)` to aggreagate the scores of each bootstrap. 

wam = function(data, num_bootstraps = 500) {
	iteration_i = 0
	
	# bootstrapping `num_bootstraps` times
	bootstraps = replicate(num_bootstraps, simplify = F, expr = { # TODO (Ahmad Alsaleh): try to use multi-core processing for replicate()
		iteration_i <<- iteration_i + 1
		single.wam_(data, iteration_i, num_bootstraps)
	}) %>% abind::abind(along = 3)
	
	# aggreagating scores of all bootstraps using the arithmatice mean
	aggregated.bootstraps = apply(bootstraps, 1:2, mean)
	
	return(aggregated.bootstraps)
}
