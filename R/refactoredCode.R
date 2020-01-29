# TODO: Add comment
# 
# Author: Rob Challen
###############################################################################

#' calculate a pointwise mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "KWindow","KNN","Kernel","Quantile"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousPointwiseConditionalEntropy = function(df, discreteVars, continuousVar, method="KWindow", ...) {
	stop("Not supported")
	switch (method,
			KWindow = calculateDiscreteContinuousPointwiseConditionalEntropy_KWindow(df, discreteVars, {{continuousVar}}, ...),
			KNN = calculateDiscreteContinuousPointwiseConditionalEntropy_KNN(df, discreteVars, {{continuousVar}}, ...),
			Kernel = calculateDiscreteContinuousPointwiseConditionalEntropy_Kernel(df, discreteVars, {{continuousVar}}, ...),
			Entropy = calculateDiscreteContinuousPointwiseConditionalEntropy_Entropy(df, discreteVars, {{continuousVar}}, ...),
			# DiscretiseByRank = calculateDiscreteContinuousMI_DiscretiseByRank(df, discreteVars, {{continuousVar}}, ...),
			# DiscretiseByValue = calculateDiscreteContinuousMI_DiscretiseByValue(df, discreteVars, {{continuousVar}}, ...),
			# Compression = calculateDiscreteContinuousMI_Compression(df, discreteVars, {{continuousVar}}, ...),
			{stop(paste0(method," not a valid option"))}
	
	)
}


#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y). I.e. the self information of X conditioned on each of the possible values for X
#' using a sliding window and local entropy measure
#' 
#' This is based on the technique described here:
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
calculateDiscreteContinuousPointwiseConditionalEntropy_KWindow = function(df, discreteVars, continuousVar, k_05=4L, ...) { 
	stop("not supported")
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
	
	# average m_i over 3 window sizes
	# k = k_05*2,
	# m_i = floor((
	#				lead(rank,n=local(k_05))-lag(rank,n=local(k_05))+
	#				lead(rank,n=local(k_05+1L))-lag(rank,n=local(k_05+1L))+
	#				lead(rank,n=local(k_05-1L))-lag(rank,n=local(k_05-1L))
	#				)/3*2)/2
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
	
	tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,!!!discreteVars) %>% summarize(
			N = max(N,na.rm = TRUE),
			N_x = max(N_x,na.rm = TRUE),
			p_x = max(as.double(N_x)/N, na.rm = TRUE),
			I_given_x = mean(I_i,na.rm = TRUE),
			I_given_x_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE))
	) %>% mutate(method = "KWindow")
	
	
	return(tmp5 %>% group_by(!!!grps))
}

#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y). I.e. the self information of X conditioned on each of the possible values for X
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousPointwiseConditionalEntropy_Kernel = function(df, discreteVars, continuousVar, collect=FALSE, ...) {
	stop("Not supported")
	df = collectDf(df,collect)
	grps = df %>% groups()
	continuousVar = ensym(continuousVar)
	
	tmp = df %>% 
			mutate(y_continuous=!!continuousVar) %>%
			group_by(!!!grps) %>%
			mutate(y_min = min(y_continuous),y_max = max(y_continuous))
	
	tmp2 = tmp %>%  
			group_by(!!!grps,!!!discreteVars) %>% 
			probabilitiesFromContinuous_Kernel(y_continuous,y_min,y_max) %>%
			rename(N_x = N, p_y_given_x=p_x, I_y_given_x = I_x)
	
	tmp3 = tmp %>%
			group_by(!!!grps) %>%
			probabilitiesFromContinuous_Kernel(y_continuous,y_min,y_max) %>%
			rename(N = N, p_y = p_x, I_y = I_x)
	
	joinList = tmp %>% joinList(defaultJoin = "y_continuous")
	
	tmp4 = tmp2 %>% 
			inner_join(tmp3, by=joinList) %>% 
			mutate(
					p_x = as.double(N_x)/N,
					pmi_xy = p_x*p_y_given_x*log(p_y_given_x/p_y)
			)
	
	tmp5 = tmp4 %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
					# integrate over dy ( well does the trapeziod part of the integration
					# TODO: this integegration could be improved. lead(pmi_xy,default = 0)
					d_I_d_xy = (pmi_xy+lag(pmi_xy,1,default=0))*as.double(y_continuous-lag(y_continuous))/2
			) %>% group_by(!!!grps,!!!discreteVars) %>% summarise (
					N = max(N,na.rm = TRUE),
					N_x = max(N_x,na.rm = TRUE),
					p_x = max(p_x,na.rm = TRUE),
					I_given_x = sum(d_I_d_xy,na.rm=TRUE)/p_x, # this does the sum over x and computes the integral over y at the same time.
					I_given_x_sd = NA,
					method = "Kernel"
			)
	
	return(tmp5)
	
}

#### ----

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a compression algorithm based on.
#' 
#' Universal and accessible entropy estimation using a compression algorithm Ram Avinery, Micha Kornreich, Roy Beck https://arxiv.org/abs/1709.10164
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_Compression = function(df, discreteVars, continuousVar, collect=FALSE, ...) {
	stop("Not supported")
	df = collectDf(df,collect)
	grps = df %>% groups()
	continuousVar = ensym(continuousVar)
	
	tmp = df %>% rename(
			y_continuous=!!continuousVar
	) 
	
	if (length(grps)==0) {
		grpsList = NULL
	} else {
		grpsList = sapply(grps,as.character)
	}
	# list of join variables for join by value
	join2List = c(grpsList, "join")
	
	# define the discretisation strategy from 4 to 256 bins depending on sample size
	# Cy defines the maximum entropy possible in that discretised 
	# map continuous variable to a stream of bytes
	
	tmp2 = tmp %>% group_by(!!!grps) %>% discretise_ByRank(y_continuous,y_discrete,binStrategy = linearBySize(8,4,256)) %>% mutate(
			y_raw = as.raw(y_discrete-1)
	)
	
	# https://www.sqlshack.com/compression-decompression-functions-sql-server-2016/
	# private function
	compressionEntropy = function(d,g,...) {
		Cy = n_distinct(d$y_discrete)
		#classes = unique(d$y_discrete)
		N = length(d$y_raw)
		C0 = length(memCompress(as.raw(rep(0,N))))
		# https://math.stackexchange.com/questions/1406747/whats-theoretical-maximum-information-compression-rate
		# in this case our we are looking to find L (compressed message length) based on max possible entropy for message type (log(Cy)) and uncompressed message length (N)
		C1 = as.double(C0+N)/log(Cy) #max(sapply(c(1:10),function(i) length(memCompress(as.raw(sample(classes,size=N,replace=TRUE)-1)))))
		C = length(memCompress(as.vector(d$y_raw)))
		if (C > C1) C=C1 # prevent entropy exceeding theoretical maximum
		H = as.double(C-C0)/(C1-C0) * log(Cy) # original paper includes a degrees of freedom parameter here. with this setup this can only be one...?
		#TODO: the log Cy term here is implicated in differences between this and other algorithms
		return(tibble(
						H=H,
						N=N
				))
	}
	
	# overall entropy of Y
	tmp3 = tmp2 %>% group_by(!!!grps) %>% group_modify(compressionEntropy) %>% rename(Hy = H) %>% mutate(join = 1)
	# Hy = calculateDiscreteEntropy(y_raw, method=method),
	# Hy = compressionEntropy(y_raw),
	# Cy = n_distinct(y_raw)),
	# Cy = max(as.integer(y_raw))-min(as.integer(y_raw))+1,
	# Cy = 256, # y is continuous therefore can have any number of classes
	# N = n()
	# ) %>% mutate(join=1)
	
	# entropy of Y given X=x
	tmp4 = tmp2 %>% group_by(!!!grps,!!!discreteVars) %>% group_modify(compressionEntropy) %>% rename(NX = N, Hygivenx = H) %>% 
			group_by(!!!grps) %>% mutate(N = sum(NX))
	
	# NX = n(),
	# Hygivenx = compressionEntropy(y_raw)
	# Hygivenx = calculateDiscreteEntropy(y_raw,Cy,method)
	# ) %>% mutate(join=1)
	
	tmp5 =  tmp4 %>% summarise(
			Hygivenx = sum( as.double(NX) / N * Hygivenx )
	) %>% mutate(join = 1)
	
	# https://en.wikipedia.org/wiki/Conditional_entropy
	# H(Y|X) = SUM(over all x) P(x)*H(Y|X=x)
	# I(X,Y) = H(Y) - H(Y|X)
	tmp6 = tmp3 %>% inner_join(tmp5, by=join2List) %>% ungroup() %>% group_by(!!!grps,N)  %>% summarise(
			I = Hy - Hygivenx,
			I_sd = as.numeric(NA),
			method = "Compression"
	)
	
	return(tmp6)
	
}

#' calculate conditional mutual information between a discrete value (X) and a continuous value (Y) using estimates of differential entropy
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @param entropyMethod - the method used to calculate the entropy (see ?tidyinfostats::calculateDiscreteEntropy) - defaults to "Grassberger"
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteContinuousPointwiseConditionalEntropy_Entropy = function(df, discreteVars, continuousVar, entropyMethod="Quantile", ...) {
	stop("Not supported")
	grps = df %>% groups()
	innerJoinList = df %>% joinList(discreteVars)
	
	continuousVar = ensym(continuousVar)
	
	tmp = df %>% 
			group_by(!!!grps, !!!discreteVars) %>% 
			calculateContinuousEntropy(!!continuousVar, method = entropyMethod, ...) %>% 
			rename(I_given_x = I, I_given_x_sd = I_sd)
	
	tmp_p_x = df %>% group_by(!!!grps) %>% 
			groupwiseCount(discreteVars, summarise = TRUE) %>% 
			mutate(p_x=as.double(N_x)/N)
	
	tmp2 = tmp %>% left_join(tmp_p_x, by=innerJoinList) %>% mutate(tmp_join = 1L)
	
	
	tmp5 = tmp2 %>% 
			group_by(!!!grps,!!!discreteVars) %>% 
			mutate(
					I_given_x = I_y*p_x - I_given_x,
					I_given_x_sd = I_y_sd*p_x + I_given_x_sd,
					method = paste0("Entropy - ",entropyMethod)
			)
	browser()
	return(tmp5)
	
}


#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
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
calculateDiscreteContinuousPointwiseConditionalEntropy_KNN = function(df, discreteVars, continuousVar, k_05=4L, useKWindow = TRUE,...) {
	grps = df %>% groups()
	maxN = df %>% group_by(!!!grps) %>% summarise(N = n()) %>% summarise(maxN = max(N, na.rm=TRUE)) %>% pull(maxN)
	if(maxN > 500 && useKWindow) {
		message("using KWindow method instead of KNN for larger sample")
		return(calculateDiscreteContinuousPointwiseConditionalEntropy_KWindow(df, discreteVars, {{continuousVar}}, k_05))
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
	
	tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,!!!discreteVars) %>% summarize(
			N = max(N,na.rm = TRUE),
			N_x = max(N_x,na.rm = TRUE),
			p_x = max(as.double(N_x)/N,na.rm = TRUE),
			I_given_x = mean(I_i,na.rm = TRUE),
			I_given_x_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE))
	) %>% mutate(method = "KNN")
	
	return(tmp5)
}




