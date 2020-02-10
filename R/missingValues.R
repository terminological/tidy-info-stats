# TODO: Add comment
# 
# Author: Rob Challen
###############################################################################

#' Adjust mutual information calculation when values are missing.
#' 
#' This corrects normal mutual information calculations for information carried by the absense of a variable. This is relevant for 
#' sparse data sets with many features such as NLP terms. Unequal "missingness" of features can contain information about the outcomes.
#' To deterimine whether a feature is missung for a given sample we need to make some assumptions. These are generally calculated
#' from the data but can alternatively be specificied directly.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (e.g. outcome)
#' @param sampleVars - the column(s) which uniquely identify the sample (e.g. person identifier)
#' @param mutualInformationFn - the function that will calculate the unadjusted MI
#' @param ... - the other parameters are passed onto the function specified in mutualInformationFn and the observedVersusExpected(...) function. Particularly sampleCount or sampleCountDf
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
adjustMIForAbsentValues = function(df, discreteVars, sampleVars, mutualInformationFn, ...) {
	grps = df %>% groups()
	outerJoinCols = df %>% joinList(discreteVars)
	innerJoinCols = df %>% joinList(defaultJoin = "tmp_join")
	
	baselineMI = mutualInformationFn(df, discreteVars, ...) %>% mutate(tmp_join = 1L) %>% compute()
	
	tmp3 = df %>% calculateDiscreteAbsentValuesMI(discreteVars, sampleVars, ...) %>%
	  rename(I_given_abs = I) %>% mutate(tmp_join=1) %>% compute()
	
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

#' calculate mutual information between a categorical value (X) and its absence in a data set.
#' 
#' This calculates the mutual information of a feature not being present in all samples
#' 
#' This is relevant for sparse data sets with many features such as NLP terms, where a term as a
#' feture may not be present in a given document, and this absense may be assymetrically distributed between
#' different classes. 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of variable (features)
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (e.g. outcome)
#' @param sampleVars - the column(s) of the sample identifier
#' @param sampleCount - an integer containing the count of all samples per outcome (discreteVars)
#' @param sampleCountDf - a dataframe containing columns for df grouping (features), and discreteVars (outcomes), N and N_x columns with expected counts of outcomes
#' see expectSamplesByOutcome(...)
#' @return a dataframe containing the distinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
observedVersusExpected = function(df, discreteVars, sampleVars, sampleCount=NULL, sampleCountDf=NULL, ...) {
  if (identical(sampleVars,NULL) && identical(sampleCount,NULL) && identical(sampleCountDf,NULL)) stop("one of sampleVars, sampleCount, or sampleCountDf must be specified")
  grps = df %>% groups()
  outerJoinCols = df %>% joinList(discreteVars)
  innerJoinCols = df %>% joinList(defaultJoin = "tmp_join")
  if (identical(sampleCountDf,NULL)) {
    if (identical(sampleCount,NULL)) {
      sampleCountDf = df %>% expectOnePerSample(discreteVars, sampleVars)
    } else {
      sampleCountDf = df %>% expectFixedSamples(discreteVars, sampleCount)
    }	
  }
  
  # How many samples per outcome (discreteVars - x) did we observe (for each feature (grps))?
  observedOutcomesPerFeature = df %>% select(!!!grps,!!!discreteVars,!!!sampleVars) %>% distinct() %>%
    group_by(!!!grps,!!!discreteVars) %>% summarise(N_x_obs = n()) 
  
  observedPerFeature = observedOutcomesPerFeature %>%
    group_by(!!!grps) %>% summarise(N_obs = sum(N_x_obs)) %>% mutate(tmp_join=1L)
  
  tmp = sampleCountDf %>% mutate(tmp_join=1L) %>%
    rename(N_exp = N, N_x_exp = N_x) %>% 
    left_join(
      observedOutcomesPerFeature, by=outerJoinCols
    ) %>% mutate(
      # fill in missing
      N_x_obs = ifelse(is.na(N_x_obs),0,N_x_obs)
    ) %>% left_join(
      observedPerFeature, by=innerJoinCols
    ) %>% mutate(
      N_obs = ifelse(is.na(N_obs),0,N_obs)
    )
  
  return(tmp %>% select(-tmp_join))
}

#' specify the number of expected observations for each class based on  a fixed number.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different features
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...) (usually the outcome)
#' @param sampleCount - the number of samples per outcome (equal numbers).
#' @return a dataframe with the expected values of all outcomes if there are a fixed number of samples
#' @import dplyr
#' @export
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
#' @return a dataframe with the expected values of all outcomes if there are a fixed number of samples per outcome
#' @import dplyr
#' @export
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
#' @return a dataframe with the expected values of all outcomes if there are identified samples and we expect every feature to be present once
#' @import dplyr
#' @export
expectOnePerSample = function(df, discreteVars, sampleVars, ...) {
	grps = df %>% groups()
	
	# the groups (features)
	left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
	
	# the unique count of samples per outcome regardless of feature
	# this will filter out outcomes which are never seen anywhere in the data
	right = df %>% ungroup() %>% select(!!!discreteVars,!!!sampleVars) %>% 
			distinct() %>% 
			group_by(!!!discreteVars) %>% 
			summarise(N_x = n()) %>% mutate(join=1)
	
	# the cross product of features, outcomes, and count of samples (per outcome)
	cross = left %>% left_join(right, by="join")
	
	# fill in the top level counts, of samples per feature
	cross = cross %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
	
	return(cross %>% select(-join))	
}

