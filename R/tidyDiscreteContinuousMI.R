#' calculate mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' This is specifically designed to supprt tidy data where there are many features, with associated values and outcomes in different columns of a dataframe or database table
#'
#' N.B. this result is the mutual information between feature value and outcome GIVEN that the feature is present. It does not account for missing values.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "KWindow","KNN","Discretise","Grassberger","Compression","Entropy","Quantile","PDF","SGolay","Kernel"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
#' @examples 
#' observations %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(outcome), value)
calculateDiscreteContinuousMI = function(df, discreteVars, continuousVar, method="KWindow", ...) {
	switch (method,
			KWindow = calculateDiscreteContinuousMI_KWindow(df, discreteVars, {{continuousVar}}, ...),
			KNN = calculateDiscreteContinuousMI_KNN(df, discreteVars, {{continuousVar}}, ...),
			#Discretisation methods
			Discretise = calculateDiscreteContinuousMI_Discretise(df, discreteVars, {{continuousVar}}, ...),
			Grassberger = calculateDiscreteContinuousMI_Discretise(df, discreteVars, {{continuousVar}}, mutualInfoMethod="Entropy", entropyMethod="Grassberger", ...),
			Compression = calculateDiscreteContinuousMI_Discretise(df, discreteVars, {{continuousVar}}, binStrategy = linearBySize(slope=8,minBins=4,maxBins=256), mutualInfoMethod="Entropy", entropyMethod="Compression", ...),
			#Continuous Entropy methods
			Entropy = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, ...),
			Quantile = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, entropyMethod="Quantile",...),
			#PDF methods
			PDF = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, entropyMethod="PDF",...),
			SGolay = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, entropyMethod="PDF", probabilityMethod = "SGolay", ...),
			Kernel = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, entropyMethod="PDF", probabilityMethod = "Kernel", ...),
			#PDF methods
			PDF2 = calculateDiscreteContinuousMI_PDF(df, discreteVars, {{continuousVar}}, ...),
			SGolay2 = calculateDiscreteContinuousMI_PDF(df, discreteVars, {{continuousVar}}, probabilityMethod = "SGolay", ...),
			Kernel2 = calculateDiscreteContinuousMI_PDF(df, discreteVars, {{continuousVar}}, probabilityMethod = "Kernel", ...),
			
			{stop(paste0(method," not a valid option"))}
	
	)
}






#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is based on the technique described here:
#' 
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' but with the important simplification of using the sliding window K elements wide rather than the k nearest neighbours. This is empirically shown to have little difference on larger datasets
#' and makes this algorithm simple to implement in dbplyr tables.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KWindow = function(df, discreteVars, continuousVar, k_05=4, ...) { 
	k_05 = as.integer(k_05)
	if (k_05<2L) k_05=2L
	grps = df %>% groups()
	continuousVar = ensym(continuousVar)
	joinList = df %>% joinList(discreteVars)
	
	df = df %>% select(!!!grps,!!!discreteVars,!!continuousVar)
	# this is confusing because groups mean 2 things here - the 
	# different types of Y (grps) which should be preserved and the categorical X 
	# has group counts (N) and subgroup counts (N_x) 
	
	tmp = df %>% groupwiseCount(discreteVars) %>% mutate(y_continuous=!!continuousVar)
	
	# the knn approach without using neighbours - i.e. a k wide sliding window
	tmp4 = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number())
	tmp4 = tmp4  %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
			
			# correct k for tails of distributions exclusive
			# kRank = row_number(),
			# m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1)+1L,
			# k = lead(kRank,n=k_05,default=max(N_x))-lag(kRank,n=k_05,default=1)+1L
			
			# correct k for tails of distributions inclusive
			# kRank = row_number(),
			# m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1),
			# k = lead(kRank,n=k_05,default=max(N_x))-lag(kRank,n=k_05,default=1)
			
			# dont correct k & exclude tails
			k = local(k_05)*2L,
			m_i = ifelse(N_x < local(k_05)*2L+1L, NA, lead(rank,n=local(k_05))-lag(rank,n=local(k_05)))
	
	)  
	
	# there are situations where this estimate can go negative - its to do with the estimation error in small sample sizez
	# and the fact that in the window method the sliding is potentially centred around different points
	# TODO: try and understand this:
	# tmp4 = tmp4 %>% mutate(m_i = ifelse(m_i < k,NA,m_i)) %>% compute()
	tmp4 = tmp4 %>% compute()
	
	# TODO: what happens when there are not enough points in one class? As is we seem to get conditional estimate based on one class
	# but not the other. In the extreme case where there are no points in the minority class we get a degnerate behaviour and
	# an overall mutual information of 0 (which may be technically correct)
	
	
	tmp4 = tmp4 %>% 
			calculateDigamma(k,digammak) %>% 
			calculateDigamma(N,digammaN) %>% 
			calculateDigamma(N_x,digammaN_x) %>% 
			calculateDigamma(m_i,digammam_i) %>% 
			mutate(
					I_i = digammaN-digammaN_x+digammak-digammam_i
			)
	
	tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,N) %>% summarize(
			I = mean(I_i,na.rm = TRUE),
			I_sd = sd(I_i,na.rm = TRUE)
	) %>% mutate(method = "KWindow")
	
	return(tmp5)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a discretisation and calculateDiscreteDiscreteMI()
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param discretiseMethod - What method will be used to discretise the data? (ByRank, ByValue, Manual)
#' @param mutualInfoMethod - What method will be used to calculate the MI once discretised?
#' @param ... - other parameters passed onto discretisation and mutual info methods
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_Discretise = function(df, discreteVars, continuousVar, discretiseMethod="ByValue", mutualInfoMethod="Histogram", ...) {
	continuousVar = ensym(continuousVar)
	grps = df %>% groups()
	
	tmp = df %>% rename(
			y_continuous=!!continuousVar
	)
	
	tmp = tmp %>% discretise(y_continuous, y_discrete, method=discretiseMethod, ...) %>% compute()
	
	return(tmp %>% calculateDiscreteDiscreteMI(discreteVars, vars(y_discrete), method=mutualInfoMethod, ...) %>% 
					collect() %>% mutate(method = paste0("Discretise - ",discretiseMethod," - ",mutualInfoMethod)))
	
}


#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is an implementation of the technique described here:
#' 
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' 
#' But it is very slow. Empirically it also does not give any better estimate that the KWindow method.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @param useKWindow - will switch to using the much faster KWindow estimator for larger sample sizes (>500) when the difference between the 2 methods is negligable
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KNN = function(df, discreteVars, continuousVar, k_05=4L, useKWindow = TRUE,...) {
	grps = df %>% groups()
	
	maxN = df %>% group_by(!!!grps) %>% summarise(N = n()) %>% summarise(maxN = max(N, na.rm=TRUE)) %>% pull(maxN)
	
	if(maxN > 500 && useKWindow) {
		message("using KWindow method instead of KNN for larger sample")
		return(calculateDiscreteContinuousMI_KWindow(df, discreteVars, {{continuousVar}}, k_05))
	}
	
	k_05 = as.integer(k_05)
	if (k_05<2L) k_05=2L
	continuousVar = ensym(continuousVar)
	
	joinList = df %>% joinList(discreteVars)
	
	df = df %>% select(!!!grps,!!!discreteVars,!!continuousVar)
	# this is confusing because groups mean 2 things here - the 
	# different types of Y (grps) which should be preserved and the categorical X 
	# NX has group counts (N) and subgroup counts (N_x) 
	tmp = df %>% groupwiseCount(discreteVars) %>% labelGroup(discreteVars, x_discrete) %>% mutate(y_continuous=!!continuousVar)
	
	# add row numbers of sequence - sorted by value, and sorted by value in groups
	tmp = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number()) %>% 
			group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(groupRank = row_number()) %>% compute()
	
	tmp = tmp %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
			rankMax = lead(rank,n=k_05*2L),
			rankMin = lag(rank,n=k_05*2L,1L)
	) %>% mutate(rankMax = ifelse(is.na(rankMax),N,rankMax))
	
	# list of join variables for join by value
	join2List = df %>% joinList(defaultJoin="join")
	
	# rhs of join by value
	tmp_join = tmp %>% ungroup() %>% rename(
					y_continuous_knn = y_continuous,
					x_discrete_knn = x_discrete,
					rank_knn = rank) %>%
			select(!!!grps,x_discrete_knn,y_continuous_knn,rank_knn) %>% mutate(join=1L) %>% compute()
	
	# TODO: this is unuseably inefficient
	
	tmp4 = tmp %>% mutate(join=1) %>% inner_join(tmp_join, by=join2List) %>% 
			filter(rankMax >= rank_knn & rankMin <= rank_knn) %>%
			mutate(y_diff = abs(y_continuous - y_continuous_knn)) %>% 
			group_by(!!!grps,rank) %>% 
			arrange(y_diff) %>% 
			mutate(sameGroup=ifelse(x_discrete_knn==x_discrete,1L,0L)) %>% #, differentGroup=ifelse(x_discrete_knn==x_discrete,0,1)) %>%
			mutate(kDist = cumsum(sameGroup), m_i = row_number()) %>%
			filter(kDist == local(k_05*2L+1L) & sameGroup==1L) %>% 
			mutate(
					k = ifelse(N < k_05*2L+1L, NA, kDist-1L), 
					m_i = m_i-1
			) %>% compute() # this is the strange definintion of knn in the original paper
	
	tmp4 = tmp4 %>% 
			calculateDigamma(k,digammak) %>% 
			calculateDigamma(N,digammaN) %>% 
			calculateDigamma(N_x,digammaN_x) %>% 
			calculateDigamma(m_i,digammam_i) %>% 
			mutate(
					I_i = digammaN-digammaN_x+digammak-digammam_i
			)
	
	tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,N) %>% summarize(
			I = mean(I_i,na.rm = TRUE),
			I_sd = sd(I_i,na.rm = TRUE)
	) %>% mutate(method = "KNN")
	
	return(tmp5)
}




#' calculate mutual information between a discrete value (X) and a continuous value (Y) using estimates of differential entropy
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @param entropyMethod - the method used to calculate the entropy (see ?tidyinfostats::calculateDiscreteEntropy) - defaults to "Grassberger"
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteContinuousMI_Entropy = function(df, discreteVars, continuousVar, entropyMethod="Quantile", ...) {
	grps = df %>% groups()
	outerJoinList = df %>% joinList(defaultJoin = "tmp_join")
	
	H_y = df %>% group_by(!!!grps) %>% 
			calculateContinuousEntropy(!!continuousVar, method = entropyMethod, ...) %>% 
			rename(I_y = I, I_y_sd = I_sd) %>% mutate(tmp_join = 1L)
	
	H_y_given_x = df %>% 
			group_by(!!!grps, !!!discreteVars) %>% 
			calculateContinuousEntropy(!!continuousVar, method = entropyMethod, ...) %>% 
			rename(N_x = N, I_given_x = I, I_given_x_sd = I_sd) %>% mutate(tmp_join = 1L)
	
	tmp2 = H_y %>% left_join(tmp_H_y, by=outerJoinList) %>% select(-tmp_join)
	
	tmp5 = tmp2 %>% 
			group_by(!!!grps, N, I_y, I_y_sd) %>% 
			summarise(
					I = I_y - sum(p_x*I_given_x),
					I_sd = I_y_sd + sum(p_x*I_given_x_sd),
					method = paste0("Entropy - ",entropyMethod)
			) %>% 
			select(-c(I_y,I_y_sd))
	
	return(tmp5)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using an estimator of PDF
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_PDF = function(df, discreteVars, continuousVar, probabilityMethod = "Kernel", ...) { 
	grps = df %>% groups()
	continuousVar = ensym(continuousVar)
	
	tmp = df %>% 
			mutate(y_continuous=!!continuousVar) %>%
			group_by(!!!grps) %>%
			mutate(y_min = min(y_continuous),y_max = max(y_continuous))
	
	tmp2 = tmp %>%  
			group_by(!!!grps,!!!discreteVars) %>% 
			probabilitiesFromContinuous(y_continuous, minVar=y_min, maxVar=y_max, method=probabilityMethod, ...) %>%
			rename(N_x = N, p_y_given_x=p_x, I_y_given_x = I_x)
	
	tmp3 = tmp %>%
			group_by(!!!grps) %>%
			probabilitiesFromContinuous(y_continuous, minVar=y_min, maxVar=y_max, method=probabilityMethod, ...) %>%
			rename(N = N, p_y = p_x, I_y = I_x)
	
	joinList = tmp %>% joinList(defaultJoin = "y_continuous")
	
	tmp4 = tmp2 %>% 
			inner_join(tmp3, by=joinList) %>% 
			mutate(
					p_x = as.double(N_x)/N,
					pmi_xy = p_x*p_y_given_x*log(p_y_given_x/p_y)
			)
	
	tmp5 = tmp4 %>% group_by(!!!grps,N,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
			# TODO: this integegration could be improved. lead(pmi_xy,default = 0)
			d_I_d_xy = (pmi_xy+lag(pmi_xy,1,default=0))*as.double(y_continuous-lag(y_continuous))/2
	) 
	
	tmp6 = tmp5 %>% group_by(!!!grps,N) %>% summarise (
			I = sum(d_I_d_xy,na.rm=TRUE), # this does the sum over x & y
			I = NA,
			method = paste0("PDF - ",probabilityMethod)
	)
	return(tmp6)
}
