
#' calculate self information when an observation has a discrete value (X).
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable - the grouping may be w.g. a test or concept. 
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param countVar - optional the column of the count variable - how often does the event happen? If missing then this will be the assumed to be individual observations. In this case the df is a contingency table
#' @param method - the method employed - valid options are "MontgomerySmith", "Histogram", "Grassberger", "InfoTheo", "Compression"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteBinaryMI =  function(df, discreteVars, countVar=NULL, method="Grassberger", ...) {
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  grps = df %>% groups()
  
  return(calculateDiscreteEntropy(df, discreteVars, countVar=!!countVar, method = method,...))
  
  # https://en.wikipedia.org/wiki/Pointwise_mutual_information
  # https://en.wikipedia.org/wiki/Information_content
  
}

#' calculate mutual information between a categorical value (X) and its absence in a data set.
#' 
#' This calculates the mutual information of a feature not being present in every sample
#' 
#' This is relevant for sparse data sets with many features such as NLP terms, where a term as a
#' feture may not be present in a given document, and this absense may be assymetrically distributed between
#' different classes. 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of variable (features)
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (e.g. outcome)
#' @param sampleVars - the column(s) of the sample identifier
#' @param sampleCount - (optional) an integer containing the count of all samples per outcome (discreteVars)
#' @param sampleCountDf - (optional) a dataframe containing columns for df grouping (features), and discreteVars (outcomes), N and N_x columns with expected counts
#' see expectSamplesByOutcome(...)
#' @return a dataframe containing the distinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteAbsentValuesMI = function(df, discreteVars, sampleVars, sampleCount=NULL, sampleCountDf=NULL, ...) {
	grps = df %>% groups()
	
	tmp = df %>% observedVersusExpected(discreteVars=discreteVars, sampleVars=sampleVars, sampleCount=sampleCount, sampleCountDf=sampleCountDf, ...)
	
	# this calculation uses the for of I(X,Y) given here:
	# https://en.wikipedia.org/w/index.php?title=Mutual_information&section=9#Relation_to_Kullback%E2%80%93Leibler_divergence
	# substituting Abs for X and our X for Y - = p_X_given_Y=y becomes p_abs_given_X=x 
	tmp2 = tmp %>% mutate(
			# we need to know I_given_unlabelled_and_x
	    p_x = as.double(N_x_exp)/ N_exp,
			p_abs = as.double(N_exp-N_obs)/N_exp,
			p_abs_given_x = as.double(N_x_exp-N_x_obs)/N_x_exp,
			N_abs = N_exp-N_obs,
			I_abs_given_x = ifelse(p_abs_given_x <= 0, 0 ,p_abs_given_x*log(p_abs_given_x/p_abs))
	)
	#tmp2: I_given_abs_and_x is pointwise mut info? or I(Abs|X=x)
	#tmp2 could be maybe combined with pointwise I(Y|X=x) to get a measure of how much
	#information Y or its absence infers about X=x 
	
	tmp3 = tmp2 %>% group_by(!!!grps,N_exp,N_abs) %>% summarise(
			# combine labelled and unlabelled I (independent)
			I = sum(p_x*I_abs_given_x),
			I_sd = as.double(NA)
	) %>% mutate(method = "Absent values") 
	
	return(tmp3)
}
  

#' calculate mutual information between a categorical value (X) and its presence in a data set.
#' 
#' This calculates the mutual information of a feature not being present in all samples
#' 
#' This is relevant for sparse data sets with many features such as NLP terms, where a term as a
#' feature is only flagged as present in the samples. 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of variable (features)
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (e.g. outcome)
#' @param sampleVars - the column(s) of the sample identifier
#' @param sampleCount - (optional) an integer containing the count of all samples per outcome (discreteVars)
#' @param sampleCountDf - (optional) a dataframe containing columns for df grouping (features), and discreteVars (outcomes), N and N_x columns with expected counts
#' see expectSamplesByOutcome(...)
#' @return a dataframe containing the distinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscretePresentValuesMI = function(df, discreteVars, sampleVars, sampleCount=NULL, sampleCountDf=NULL, ...) {
	
  grps = df %>% groups()
	tmp = df %>% observedVersusExpected(discreteVars=discreteVars, sampleVars=sampleVars, sampleCount=sampleCount, sampleCountDf=sampleCountDf, ...)
	
	# this calculation uses the for of I(X,Y) given here:
	# https://en.wikipedia.org/w/index.php?title=Mutual_information&section=9#Relation_to_Kullback%E2%80%93Leibler_divergence
	# substituting Pres for X and our X for Y - = p_X_given_Y=y becomes p_abs_given_X=x 
	tmp2 = tmp %>% mutate(
			# we need to know I_given_unlabelled_and_x
	    p_x = as.double(N_x_exp)/ N_exp,
			p_pres = as.double(N_obs)/N_exp,
			p_pres_given_x = as.double(N_x_obs)/N_x_exp,
			I_pres_given_x = ifelse(p_pres_given_x <= 0, 0 ,p_pres_given_x*log(p_pres_given_x/p_pres))
	)
	#tmp2: I_given_abs_and_x is pointwise mut info? or I(Abs|X=x)
	#tmp2 could be maybe combined with pointwise I(Y|X=x) to get a measure of how much
	#information Y or its absence infers about X=x 
	
	tmp3 = tmp2 %>% group_by(!!!grps,N_exp,N_obs) %>% summarise(
			# combine labelled and unlabelled I (independent)
			I = sum(p_x*I_pres_given_x),
			I_sd = as.double(NA)
	) %>% mutate(method = "Present values") 
	
	return(tmp3)
}