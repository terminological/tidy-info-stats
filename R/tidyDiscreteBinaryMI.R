
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
  
  return(calculateEntropy(df, discreteVars, countVar=!!countVar, method = method,...))
  
  # https://en.wikipedia.org/wiki/Pointwise_mutual_information
  # https://en.wikipedia.org/wiki/Information_content
  
}


#' calculate self information when an observation has a discrete value (X).
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable - the grouping may be w.g. a test or concept. 
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param countVar - optional the column of the count variable - how often does the event happen? If missing then this will be the assumed to be individual observations. In this case the df is a contingency table
#' @param method - the method employed - valid options are "Histogram", "Grassberger"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateSelfInformation =  function(df, discreteVars, countVar=NULL, method="Histogram", ...) {
  switch (method,
          Histogram = calculateSelfInformation_Histogram(df, groupVars, countVar={{countVar}}, ...),
          Grassberger = calculateSelfInformation_Grassberger(df, groupVars, countVar={{countVar}}, ...),
          {stop(paste0(method," not a valid option"))}
  )
}
  
  
  
