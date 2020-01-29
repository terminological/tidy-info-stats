# TODO: Add comment
# 
# Author: Rob Challen
###############################################################################

#' estimate PDFs from samples in a tidy friendly manner
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param continuousVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "Kernel", "SGolay"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
probabilitiesFromContinuous = function(df, continuousVar, method="Kernel", ...) {
	switch (method,
			Kernel = probabilitiesFromContinuous_Kernel(df, {{continuousVar}}, ...),
			SGolay = probabilitiesFromContinuous_SGolay(df, {{continuousVar}}, ...),
			{stop(paste0(method," not a valid option"))}
	
	)
}


#' Helper function to calculate probability from continuous data in a tidy friendly manner
#'
#' The purpose of this is to calculate the probabilities of events from continuous data. 
#' This function is useful when you have a set of observations from a continuous distribution.
#' 
#' @param df a dataframe containing a column of a continuous variable X and one row per observation, 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param continousVar the datatable column(s) containing the observation.
#' @param k_05 the sgolay smoothing window width.
#' @return A mutated datatable with observations of X, the total number of observations of X (N), the probability density (p_x)
#' @import dplyr
#' @export
probabilitiesFromContinuous_SGolay = function(df, continuousVar, k_05 = 10) {
	continuousVar = ensym(continuousVar)
	grps = df %>% groups()
	
	# groupwise count creates an N and N_x  column based on groupVars, and countVar
	df = df %>% mutate(tmp_x_continuous = !!continuousVar)
	tmp2 = df %>%  group_by(!!!grps) %>% applySGolayFilter(tmp_x_continuous, d_x_d_r, k_05 = k_05, p=2, m=1)
	tmp2 = tmp2 %>% mutate(
	    p_x = 1.0/ifelse(d_x_d_r <= 0, NA, d_x_d_r)
	) %>% select( -d_x_d_r )
	
	tmp2 = tmp2 %>% rename(!!continuousVar := tmp_x_continuous)
	
	return(tmp2)
}

#' Helper function to calculate probability from continuous data in a tidy friendly manner
#'
#' The purpose of this is to calculate the probabilities of events from continuous data. 
#' This function is useful when you have a set of observations from a continuous distribution.
#' Kernel methods resample the data and produce an evenly spaced output function so needs some sort of support range defined for each group (minVar and maxVar)
#' 
#' @param df a dataframe containing a column of a continuous variable X and one row per observation, 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param continousVar - the datatable column(s) containing the observation.
#' @param minVar - the name of the column containing the minimum value for the observation / support range
#' @param maxVar - the name of the column containing the maximum value for the observation / support range
#' @return A mutated datatable with observations of X, the total number of observations of X (N), the probability density (p_x), and self information (I_x) associated with the value of X
#' @import dplyr
#' @export
probabilitiesFromContinuous_Kernel = function(df, continuousVar, minVar=NULL, maxVar=NULL, collect=FALSE) {
	df = collectDf(df,collect)
	continuousVar = ensym(continuousVar)
	minVar = tryCatch(ensym(minVar),error = function(e) NULL)
	maxVar = tryCatch(ensym(maxVar),error = function(e) NULL)
	grps = df %>% groups()
	
	# groupwise count creates an N and N_x  column based on groupVars, and countVar
	df = df %>% mutate(tmp_x_continuous = !!continuousVar, N = n())
	
	if (!identical(minVar,NULL)) {
		df = df %>% mutate(tmp_x_min = !!minVar)
	} else {
		df = df %>% mutate(tmp_x_min = min(tmp_x_continuous))
	}
	
	if (!identical(maxVar,NULL)) {
		df = df %>% mutate(tmp_x_max = !!maxVar)
	} else {
		df = df %>% mutate(tmp_x_max = max(tmp_x_continuous, na.rm = TRUE))
	}
	
	tmp2 = df %>% group_by(!!!grps) %>% group_modify(
			function(d,...) {
				dens = density(d$tmp_x_continuous, bw="nrd0", kernel="gaussian", from=min(d$tmp_x_min,na.rm = TRUE), to=max(d$tmp_x_max,na.rm = TRUE), n=512)
				originalSize = max(d$N,na.rm = TRUE)
				return(
						tibble(
								N = originalSize,
								tmp_x_continuous = dens$x,
								p_x = dens$y
						)
				)
			}
	)
	
	tmp2 = tmp2 %>% rename(!!continuousVar := tmp_x_continuous)
	if (!identical(minVar,NULL)) {
		tmp2 = tmp2 %>% rename(!!minVar := tmp_x_min)
	} else {
		tmp2 = tmp2 %>% select(-tmp_x_min)
	}
	
	if (!identical(maxVar,NULL)) {
		tmp2 = tmp2 %>% rename(!!maxVar := tmp_x_max)
	} else {
		tmp2 = tmp2 %>% select(-tmp_x_max)
	}
	
	return(tmp2)
}
