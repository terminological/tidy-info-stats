#' calculate entropy of a sequence of continuous variables 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param method - the method employed - valid options are "MontgomerySmith", "Histogram", "Grassberger", "InfoTheo", "Compression"
#' @param ... - the other parameters are passed onto the implementations
#' @return a single value for the entropy of the vector
#' @import dplyr
#' @export
calculateContinuousEntropy = function(df, continuousVar, method, ...) {
  switch (method,
          Quantile = calculateContinuousEntropy_Quantile(df, {{continuousVar}}, ...),
		      PDF = calculateContinuousEntropy_PDF(df, {{continuousVar}}, ...),
          {stop(paste0(method," not a valid option"))}
  )
}

# TODO:
# If there a method involving horizontal visibility?
# https://stats.stackexchange.com/questions/97676/kozachenko-leonenko-entropy-estimation - 
# a join on dense rank + 1, diff the values, diff the not dense rank, and then apply to digamma


#' calculate differential entropy of a continuous value (X) using a quantile function smoothing approach:
#' Entropy of continuous data - e.g. SGolay approach / discretisation
#' 
#' The Savitsky Golay filter width required to make this at all accurate needs a reasonable amount of data ~ 20 points
#' S. M. Sunoj and P. G. Sankaran, “Quantile based entropy function,” Stat. Probab. Lett., vol. 82, no. 6, pp. 1049–1053, Jun. 2012, doi: 10.1016/j.spl.2012.02.005. [Online]. Available: http://www.sciencedirect.com/science/article/pii/S0167715212000521
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - the half window width of the SG filter that smooths the data. This is dependent on data but typically not less that 10.
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateContinuousEntropy_Quantile = function(df, continuousVar, k_05=10, ...) {
  grps = df %>% groups()
  
  continuousVar = ensym(continuousVar)
  #joinList = joinList(df, discreteVars)
  
  tmp = df %>% group_by(!!!grps) %>% mutate(
    N = n(),
    y_continuous = !!continuousVar
  )
  
  tmp2 = tmp %>%  group_by(!!!grps) %>% applySGolayFilter(y_continuous, d_Q_d_p, k_05 = k_05) %>% mutate(
    log_d_Q = log(ifelse(d_Q_d_p <= 0, NA, d_Q_d_p))
  ) 
  
  tmp4 = tmp2 %>% summarise(
      N = max(N,na.rm = TRUE),
      I = sum(log_d_Q*(1.0/(N+1L)),na.rm = TRUE),
      I_sd = as.double(NA),
      method="Quantile"
  )
  
  return(tmp4)
  
  
}


#' calculate differential entropy of a continuous value (X) using a probability function approach:
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param continuousVar - the column of the continuous value (Y)
#' @param probabilityMethod - the method of calculating a PDF - valid options are "SGolay", "Kernel"
#' @param ... - passed onto probabilitiesFromContinuous
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateContinuousEntropy_PDF = function(df, continuousVar, probabilityMethod="Kernel", ...) {
	
	grps = df %>% groups()
	continuousVar = ensym(continuousVar)
	
	tmp = df %>% 
			mutate(x_continuous=!!continuousVar) %>%
			group_by(!!!grps)
	
	tmp2 = tmp %>%
			probabilitiesFromContinuous(x_continuous, method=probabilityMethod, ...) %>%
			mutate(p_x_I_x = ifelse(p_x <= 0, 0, -p_x*log(p_x)))
	
	tmp3 = tmp2 %>% group_by(!!!grps, N) %>% arrange(x_continuous) %>% mutate(
			d_I_d_x = (p_x_I_x+lag(p_x_I_x,1,default=0))*as.double(x_continuous-lag(x_continuous))/2
		) %>% summarise(
			I = sum(d_I_d_x,na.rm=TRUE),
			I_sd = NA,
			method = paste0("PDF - ",probabilityMethod)
		)
	
	return(tmp3)
	
	
}