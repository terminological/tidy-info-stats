# TODO: Add comment
# 
# Author: Rob Challen
###############################################################################

#' calculate mutual information between a categorical value (X) and a continuous value (Y) where a total number of possible outcomes is known.
#' 
#' This corrects normal mutual information calculations for information carried by the absense of a variable. This is relevant for 
#' sparse data sets with many features such as NLP terms.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (e.g. outcome
#' @param mutualInformationFn - the function that will calculate the unadjusted MI
#' @param sampleCountDf - a dataframe containing colums for grouping, and discreteVars, N and N_x columns with counts
#' see expectFixedSamples(...)
#' @param ... - the other parameters are passed onto the function specified in mutualInformationFn 
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
adjustMIForAbsentValues = function(df, discreteVars, mutualInformationFn, sampleCountDf, ...) {
	grps = df %>% groups()
	outerJoinCols = df %>% joinList(discreteVars)
	innerJoinCols = df %>% joinList(defaultJoin = "tmp_join")
	
	baselineMI = mutualInformationFn(df, discreteVars, ...) %>% mutate(tmp_join = 1L)
	
	observedCount = df %>% groupwiseCount(discreteVars,summarise=TRUE) %>% rename(N_obs=N, N_x_obs=N_x)
	
	tmp3 = df %>% calculateDiscreteAbsentValuesMI(discreteVars, sampleCountDf = sampleCountDf, ...) %>%
	  rename(I_given_abs = I) %>% mutate(tmp_join=1)
	
	tmp4 = tmp3 %>% left_join(
			baselineMI %>% select(-N) %>% rename(I_given_pres = I, method_sub = method), by=innerJoinCols
		) %>% mutate(
			I_given_pres = ifelse(is.na(I_given_pres),0,I_given_pres),
			method = paste0("Absent + ",method_sub),
			I = I_given_pres+I_given_abs,
			I_sd = as.double(NA)
		)
	
	return(tmp4)
}


#' summary info about the observations
#' 
#' @param df - may be grouped, in which case the value is interpreted as different features
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param expected - the total number of samples for each outcome expected if there were no missing values.
summariseObservations = function(df, discreteVars, countVar = NULL) {
	countVar = tryCatch(ensym(countVar),error = function(e) NULL)
	grps = df %>% groups()
	outerJoin = df %>% joinList(discreteVars)
	
	left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
	right = df %>% ungroup() %>% select(!!!discreteVars) %>% distinct() %>% mutate(join=1)
	cross = left %>% left_join(right, by="join")
	
	cross = cross %>% left_join(df %>% groupwiseCount(discreteVars, !!countVar, summarise=TRUE), by=outerJoin) %>%
			mutate(
					N_x = ifelse(is.na(N_x), 0, N_x),
					N = ifelse(is.na(N), 0, N)
			)
	return(cross %>% select(-join))
}

#' specify the number of expected observations for each class based on  a fixed number.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different features
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param sampleCount - the number of samples per outcome (equal numbers).
expectFixedSamples = function(df, discreteVars, sampleCount, ...) {
	grps = df %>% groups()
	left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
	right = df %>% ungroup() %>% select(!!!discreteVars) %>% distinct() %>% mutate(join=1)
	cross = left %>% left_join(right, by="join")
	cross = cross %>% mutate(N_x = sampleCount) %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
	return(cross %>% select(-join))	
}

#' specify the number of expected observations for each class based on  a fixed number.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (usually the outcome)
#' @param samplePerOutcomeDf - a dataframe containing one row per outcome and the number of samples expected for that outcome in column N
expectSamplesByOutcome = function(df, samplesPerOutcomeDf, ...) {
	grps = df %>% groups()
	left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
	cross = left %>% left_join(samplesPerOutcomeDf %>% mutate(join=1), by="join",copy=TRUE)
	cross = cross %>% mutate(N_x = N) %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
	return(cross %>% select(-join))	
}

#' specify the number of expected observations for each class based on a grouping of the data
#' 
#' @param df - may be grouped, in which case the value is interpreted as different features
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) i.e. typically the outcome
#' @param sampleVars - the column(s) that identify the sample space, quoted by vars(...) i.e. usually the sample identifier  
expectOnePerSample = function(df, discreteVars, sampleVars, ...) {
	grps = df %>% groups()
	
	# the groups (features)
	left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
	
	# the unique count of samples per outcome regardless of feature
	# this will filter out outcomes which are never seen anywhere in the data
	right = df %>% select(!!!discreteVars,!!!sampleVars) %>% 
			distinct() %>% 
			group_by(!!!discreteVars) %>% 
			summarise(N_x = n()) %>% mutate(join=1)
	
	# the cross product of features, outcomes, and count of samples (per outcome)
	cross = left %>% left_join(right, by="join")
	
	# fill in the top level counts, of samples per feature
	cross = cross %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
	
	return(cross %>% select(-join))	
}

