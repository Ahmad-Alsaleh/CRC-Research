"
This code contains functions that will be used by a different file to compute
the AUC values of different ML models using differnt subsets of features. 

Functions ending with an underscore are private functions meant to be used in
this file only.

"
library(ggplot2)

# Returns a transformation object with reversed logarithmically spaced values ----
# Used for the x-axis of the plot returned by generatePlot()
reverselog_trans_ = function() {
	trans = function(x) -log(x, 10)
	inv = function(x) 10^(-x)
	scales::trans_new(paste0("reverselog-", format(10)), trans, inv,
										scales::log_breaks(base = 10), domain = c(-Inf, Inf))
}	

# Plots the AUC values of an ML model using ggplot2 ----
# Prerequisites: the column that contains the size of the subsets is called "subset.size"
# Arguments:
# `AUC.values`: a data.frame returned from compute.AUCs()
# `model.name`: a string contains the name of the model that was used to obtain `AUC.values`
generatePlot = function(AUC.values, model.name) {
	melted.AUC.values = reshape2::melt(AUC.values, id = "subset.size")
	plot = ggplot(melted.AUC.values, aes(x = subset.size, y = value, color = variable)) +
		xlab("Features Used") +
		ylab("AUC") +
		ggtitle(model.name) +
		theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
		scale_x_continuous(trans = reverselog_trans_()) +
		geom_line(linetype = "dotdash", linewidth = 1) +
		labs(linetype = "FS Method", color = "FS Method") + scale_y_continuous(
			breaks = seq(min(plot$data$value), max(plot$data$value), length.out = 5),
			labels = scales::label_number(accuracy = 0.01))
	
	return(plot)
}

# Computes a single AUC value for a given subset of features with a given subset size ----
compute.AUC_ = function(features, data.train, data.test, predict.func, subset.size, model.func, ...) {
	# creating formula
	formula = features[1:subset.size, ] %>% paste0(collapse = "+") %>% paste0("y~", .) %>% as.formula
	# fitting model
	model.fit = model.func(formula = formula, data = data.train, ...)
	
	# computing AUC value
	AUC = predict.func(model.fit, data.test) %>%
		ROCR::prediction(labels = data.test$y) %>%
		ROCR::performance(measure = "auc")
	AUC = AUC@y.values[[1]]
	
	return(AUC)
}

# Computes all AUC values for each set of features with different subset sizes ----

# Arguments:
# `data`: the dataset used to fit the model. Last column contains the response
#				and should be called `y`. 
# `imp.features`: a data.frame where each column contains the names of the
#				features to be used in fitting the model.
# `predict.func`: a function with two arguments: `model.fit` and `data.train` in
#				this order and returns the predictions probabilities of levels(data$y)[2].
#				If the model's predict function takes different arguments, use a wrapper
#				function. 
# `model.func`: a function with two arguments of the following names: `formula`
#				and `data`. If the model's implementation takes different arguments,
#				use a wrapper function.
# `...`: optional arguments passed to `model.func`.
# `folds`: a list of vectors where each vector contains the indices of testing samples.
#				If a single value is passed then it corresponds to the number of folds.
# `subset.sizes`: a vector containing the different number of features to be
#				used when fitting the model. If a single value is passed then it
#				corresponds to the number of subset sizes to use.

compute.AUCs = function(data, imp.features, predict.func, model.func, ..., folds = 5, subset.sizes = 10) {
	if (length(folds) == 1)
		folds = caret::createFolds(data$y, k = folds)
	
	if (length(subset.sizes) == 1)
		subset.sizes = floor(exp(seq(log(5), log(nrow(imp.features)),
																 length.out = subset.sizes)))
	
	# `AUCs` stores the AUC values of a single fold
	# rows are number of features used (the size of the subset)
	# and columns are the features selection methods
	AUCs = matrix(nrow = length(subset.sizes), ncol = ncol(imp.features),
								dimnames = list(subset.sizes, colnames(imp.features)))
	
	# setting up parallel processing for CV (Windows not supported)
	if (Sys.info()["sysname"] == "Windows") {
		cores_n = 1
		message("Multi-core processing on is not supported for Windows.\nRunning on a single core\n")
	}
	else {
		cores_n = parallelly::availableCores()
		message("Multi-core processing using ", cores_n, " core(s)\n")
	}
	
	# cross validation to find AUC values
	cv.AUCs = parallel::mclapply(folds, mc.cores = cores_n, FUN = function(fold) {
		for(feature_i in 1:ncol(imp.features)) {
			for(subset.size_i in 1:length(subset.sizes)) {
				set.seed(42)
				AUCs[subset.size_i, feature_i] =
					compute.AUC_(imp.features[feature_i], data[-fold, ], data[fold, ],
											 predict.func, subset.sizes[subset.size_i], model.func, ...)
			} # next subset of features
		} # next feature selection method
		return(AUCs %>% as.data.frame)
	})
	
	# combining AUCs values of the k folds (using the arithmetic mean)
	cv.AUCs = cv.AUCs %>% abind::abind(along = 3) %>%
		apply(1:2, mean) %>% as.data.frame
	cv.AUCs = cbind(subset.size = subset.sizes, cv.AUCs)
	rownames(cv.AUCs) = NULL
	return(cv.AUCs)
}
