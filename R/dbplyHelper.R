#' Some databases don't support window functions over grouped data. This requires a workaround to group and summarise then join the data
#' 
#' @param df - a df which may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param ... - the contents of the mutuate function
#' @return a dbplyr dataframe containing the grouped function
#' @export
group_mutate <- function(df, ...) {
  grps = df %>% groups()
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList, "tmp_join")
  
  tmp_df <- df %>% group_by(!!!grps) %>% summarize(...) %>% ungroup() %>% mutate(tmp_join=1)
  
  return(df %>% mutate(tmp_join=1) %>% left_join(tmp_df, by=joinList) %>% select(-tmp_join))
}


#' Fail if table is a dbplyr table and the user does not request collection
#' @param df - a dataframe
#' @param collect - boolean, should the table be collected if it si a dbplyr table.
#' @export
collectDf = function(df, collect) {
  if ("tbl_sql" %in% class(df)) {
    if (collect) {
      df %>% collect()
    } else {
      stop("This implementation does not support dbplyr data frames")
    }
  }
  df
}

#' a set of the first 50 digamma values
digammaLookup = tibble(n = c(1:50), digamma = digamma(c(1:50)))

#' Calculate an estimate of the digamma function value for dplyr and dbplyr dataframes. This will be done in SQL if possible.
#' 
#' @param df - a dataframe
#' @param inputVar - the reference to the column for which the digamma value will be calculated (this must be an integer)
#' @param outputVar - the column that the digamma result will be written to
#' @export
calculateDigamma = function(df, inputVar, outputVar) {
  inputVar = ensym(inputVar)
  outputVar = ensym(outputVar)
  tmp = df %>% 
    left_join(digammaLookup %>% rename(!!inputVar := n, !!outputVar := digamma), by=as.character(inputVar), copy=TRUE) %>%
    mutate(!!outputVar := ifelse(is.na(!!outputVar), log(!!inputVar-1) + 1.0/(2*(!!inputVar-1) - 1.0/(12*(!!inputVar-1)^2) ), !!outputVar))
}