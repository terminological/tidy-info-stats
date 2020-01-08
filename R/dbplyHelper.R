#' Some databases don't support window functions over grouped data. This requires a workaround to group and summarise then join the data
#' 
#' @param df - a df which may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param ... - the contents of the mutuate function
#' @return a dbplyr dataframe containing the grouped function
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