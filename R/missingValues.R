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
#' see expectSamplesByOutcome(...)
#' @param ... - the other parameters are passed onto the function specified in mutualInformationFn 
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
adjustDiscreteMIForAbsentValues = function(df, discreteVars, mutualInformationFn, sampleCountDf, ...) {
	grps = df %>% groups()
	outerJoinCols = df %>% joinList(discreteVars)
	innerJoinCols = df %>% joinList(defaultJoin = "tmp_join")
	
	baselineMI = mutualInformationFn(df, discreteVars, ...) %>% mutate(tmp_join = 1L)
	
	observedCount = df %>% groupwiseCount(discreteVars,summarise=TRUE) %>% rename(N_obs=N, N_x_obs=N_x)
	
	tmp = sampleCountDf %>% mutate(tmp_join=1L) %>%
			rename(N_exp = N, N_x_exp = N_x) %>% 
			left_join(
				observedCount %>% select(-N_obs), by=outerJoinCols
			) %>% mutate(
				N_x_obs = ifelse(is.na(N_x_obs),0,N_x_obs)
			) %>% left_join(
				# this is needed to fill N_obs for missing columns
				observedCount %>% select(!!!grps,N_obs) %>% distinct() %>% mutate(tmp_join=1L), by=innerJoinCols
			) %>% mutate(
				N_obs = ifelse(is.na(N_obs),0,N_obs)
			)
	
	
	# this calculation uses the for of I(X,Y) given here:
	# https://en.wikipedia.org/w/index.php?title=Mutual_information&action=edit&section=9
	# substituting Abs for X and our X for Y - = p_X_given_Y=y becomes p_abs_given_X=x 
	tmp2 = tmp %>% mutate(
					# we need to know I_given_unlabelled_and_x
					p_abs = as.double(N_exp-N_obs)/N_exp,
					p_abs_given_x = as.double(N_x_exp-N_x_obs)/N_x_exp,
					p_pres = 1-p_abs,
					I_abs_given_x = ifelse(p_abs_given_x <= 0, 0 ,p_abs_given_x*log(p_abs_given_x/p_abs))
			)
	#tmp2: I_given_abs_and_x is pointwise mut info? or I(Abs|X=x)
	#tmp2 could be maybe combined with pointwise I(Y|X=x) to get a measure of how much
	#information Y or its absence infers about X=x 
			
	tmp3 = tmp2 %>% group_by(!!!grps,N_exp,N_obs) %>% summarise(
					# combine labelled and unlabelled I (independent)
					I_given_abs = sum(p_abs*I_abs_given_x)
			) %>% mutate(tmp_join = 1L) 
	
	tmp4 = tmp3 %>% left_join(
			baselineMI %>% select(-N) %>% rename(I_given_pres = I), by=innerJoinCols
		) %>% mutate(
			I_given_pres = ifelse(is.na(I_given_pres),0,I_given_pres)
		)
	
	tmp5 = tmp4 %>% group_by(!!!grps, method, N_exp, N_obs) %>% summarise(
					I = I_given_pres+I_given_abs,
					# I = p_pres*I_given_pres + p_abs*I_given_abs,
					# pretty convinced probabilities not required here 
					# as they are part of the I calculations
					I_sd = as.double(NA)
			) %>% ungroup() %>% mutate(
					method = paste0("Absent + ",method)
			)
	# browser()
	return(tmp5)
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

