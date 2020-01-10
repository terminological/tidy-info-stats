
#' calculate mutual information between a discrete value (X) with a binary outcome.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable - the grouping may be w.g. a test or concept. 
#' @param discreteVar - the column of the discrete value (X)
#' @param countVar - optional the column of the count variable - how often does the event happen? If missing then this will be the assumed to be individual observations. In this case the df is a contingency table
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteBinaryMI =  function(df, discreteVar, countVar=NA, ...) {
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  grps = df %>% groups()
  discreteVar = ensym(discreteVar)
  
  # count = number of occurrences of each class
  # documents = number of documents that each concept appears in (useful for idf)
  if (identical(countVar,NULL)) {
    tmp = df %>% group_by(!!!grps, !!discreteVar) %>% summarise(count = n())
  } else {
    tmp = df %>% group_by(!!!grps, !!discreteVar) %>% summarise(count = sum(!!countVar))
  }
  
  # TODO: Not really sure where this si going.
  # https://en.wikipedia.org/wiki/Pointwise_mutual_information
  # https://en.wikipedia.org/wiki/Information_content
  
}