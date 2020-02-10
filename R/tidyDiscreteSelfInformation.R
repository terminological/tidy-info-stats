#' calculate self information of a discrete variable 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of discrete variable
#' @param discreteVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param method - the method employed - valid options are "Histogram"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the distinct values of the discrete variable, per group, with columns I_x, and method
#' @import dplyr
#' @export
calculateDiscreteSelfInformation = function(df, discreteVars, method="Histogram", ...) {
	switch (method,
		Histogram = calculateDiscreteSelfInformation_Histogram(df, discreteVars, ...),
		Grassberger = calculateDiscreteSelfInformation_Grassberger(df, discreteVars, ...),
		{stop(paste0(method," not a valid option"))}
	)
}

#' calculate the self information for all values of a discrete variable in a dplyr friendly way
#' 
#' see: R. de Matos Simoes and F. Emmert-Streib, “Influence of statistical estimators of mutual information and data heterogeneity on the inference of gene regulatory networks,” PLoS One, vol. 6, no. 12, p. e29279, Dec. 2011 [Online]. Available: http://dx.doi.org/10.1371/journal.pone.0029279
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param countVar - (optional) if this datafram represents summary counts, the columns of the summary variable.
#' @param mm - Apply a miller-madow adjustment to the result? default = TRUE
#' @return a dataframe containing the distinct values of the discrete variable, per group, with columns I_x, and method
#' @import dplyr
#' @export
calculateSelfInformation_Histogram = function(df, groupVars, countVar=NULL, mm=TRUE, ...) {
  grps = df %>% groups()
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  # groupVars = ensyms(groupVars)
  
  tmp = df %>% groupwiseCount(groupVars, !!countVar, summarise=TRUE) %>% 
    group_by(!!!grps) %>% 
    mutate(
      p_x = as.double(N_x)/N
    )
  if (mm) {
    tmp = tmp %>% mutate(
      C_x = n(),
      I_x = -log(p_x) + as.double(C_x-1)/(2*N*p_x*C_x), #the miller-madow adjustment
      method = "Histogram MM"
    ) %>% select(-C_x)
  } else {
    tmp = tmp %>% mutate(
      I_x = -log(p_x),
      method = "Histogram"
    )
  }
  
  return(tmp)
}

#' calculate self information of a discrete value (X) using a histogram approach using the following method
#' 
#' P. Grassberger, “Entropy Estimates from Insufficient Samplings,” arXiv [physics.data-an], 29-Jul-2003 [Online]. Available: http://arxiv.org/abs/physics/0307138
#' 
#' but with a digamma based function (rather than harmonics) detailed in eqns 31 & 35.
#' For our purposes we fix l=0 to give the form in eqn 27. The error in this method is supposedly better for undersampled cases (where number of bins similar to number of samples)
#' 
#' This is a bit of a cheat as works out the overall entropy and then scales that to get the self information but seems to produce the right answer
#'  
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param countVar - (optional) if this datafram represents summary counts, the columns of the summary variable.
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateSelfInformation_Grassberger = function(df, groupVars, countVar=NULL, ...) {
	grps = df %>% groups()
	countVar = tryCatch(ensym(countVar),error = function(e) NULL)
	# groupVars = ensyms(groupVars)
	
	tmp = df %>% groupwiseCount(groupVars, !!countVar, summarise=TRUE) %>% group_by(!!!grps) %>% mutate(C_x = n())
	tmp2 = tmp %>% calculateDigamma(N_x,digamma_N_x) %>% 
			mutate(
					p_x = as.double(N_x)/N,
					G_N_x = digamma_N_x + ((-1)^N_x)/(N_x*(N_x+1)), #TODO: consider expanding for l=1, l=2, etc...
					I2_x = -log(p_x)
			)
	
	tmp3 = tmp2 %>% ungroup() %>% group_by(!!!grps) %>% mutate(
			H2 = sum(p_x*I2_x),
			H = log(N) - sum(p_x * G_N_x,na.rm = TRUE),
			I_x = H/H2*I2_x
	) %>% select(-c(H2,H,I2_x,G_N_x))
	
	return(tmp3 %>% ungroup() %>% mutate(method="Grassberger"))
}
