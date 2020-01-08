#' Helper function to calculate probability from grouped data in a tidy friendly manner
#'
#' The purpose of this is to calculate the probabilities of all binary outcomes from class data. 
#' This function is useful when you have 2 types of events (X and Y) and
#' you either have a set of observations of their co-occurrence, containing non-unique X & Y combinations, or you have a confusion matrix of the counts of their combinations where 
#' each row has a unique combination of X and Y and a third column contains the counts of XY co-occurrences.
#' 
#' https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-461
#' 
#' $$I = H_x + H_y - H_xy$$
#' 
#' Where I is information
#' 
#' x is a member of X and |X| is number of possible values of x with no non empty bins (classes_X):
#' y is a member of Y and |Y| is number of possible values of y with no non empty bins (classes_Y):
#' xy is a member of XY and |XY| is number of possible values of combination of x and y with no non empty bins (classes_XY):
#' N is the number of samples
#' 
#' $$H_x__mm = H_x__emp + (|Z_x| - 1)/2N$$
#' 
#' $$I__mm = I__emp + (|X| - 1)/2N + (|Y| - 1)/2N - (|XY| - 1)/2N$$
#' 
#' $$I__mm = I__emp + (|X| + |Y| - |XY| - 1)/2N$$
#' 
#' This is a Miller-Maddow adjustment.
#' 
#' @param df a dataframe containing 2 columns defining class of event X and class of event Y and either one row per event, 
#' or a count of observations, for each class combination. 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param groupXVar the datatable column defining the class of event X
#' @param groupXVar the datatable column defining the class of event Y
#' @param countVar the datatable column containing the observed frequency combination of event XY. If this is missing the row count will be used instead
#' @return A new datatable with all possible combinations of X&Y and the probabilities associated with each outcome (i.e. an N(X) by N(Y) set of binary confusion matrices)
#' @import dplyr
#' @export
probabilitiesFromGroups = function(df, groupXVar, groupYVar, countVar=NULL) {
  grps = df %>% groups()
  groupXVar = ensym(groupXVar)
  groupYVar = ensym(groupYVar)
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList,"join")
  if (identical(countVar,NULL)) {
    # there is no count column.  The data is grouped and has event X and event Y entried for each occurrence
    df = df %>% group_by(!!!grps, !!groupXVar, !!groupYVar) %>% summarise(f_xy = n())
  } else {
    # there is a count column.  The data is grouped and has event X and event Y and counts of occurrences
    df = df %>% group_by(!!!grps, !!groupXVar, !!groupYVar) %>% summarise(f_xy = sum(!!countVar))
  }
  N = df %>% ungroup() %>% group_by(!!!grps) %>% summarise(
    N=sum(f_xy), 
    classes_XY=n_distinct(!!groupXVar,!!groupYVar),
    classes_X=n_distinct(!!groupXVar), 
    classes_Y=n_distinct(!!groupXVar)
    ) %>% 
    mutate(mm_adjust = as.double(classes_X+classes_Y-classes_XY-1)/(2*N)) %>% 
    mutate(join=1)
  # N = N %>% mutate(classes_XY = nrow(N))
  X = df %>% ungroup() %>% group_by(!!!grps,!!groupXVar) %>% summarise(f_x = sum(f_xy)) %>% mutate(join=1) # grouping
  # X = X %>% mutate(classes_X = nrow(X))
  Y = df %>% ungroup() %>% group_by(!!!grps,!!groupYVar) %>% summarise(f_y = sum(f_xy)) %>% mutate(join=1) # grouping
  #Y = Y %>% mutate(classes_Y = nrow(Y))
  # I__mm = I__emp + (|X| + |Y| - |XY| - 1)/2N
  
  
  XY = (X %>% inner_join(Y, by=joinList) %>% inner_join(N, by=joinList)) %>% select(-join)
  joinAll = c(grpsList,as.character(groupXVar),as.character(groupYVar))
  XY = XY %>% 
    left_join(df, by=joinAll) %>% 
    mutate( f_xy = ifelse(is.na(f_xy),0,f_xy))
  return( XY %>% group_by(!!!grps) %>% probabilitiesFromCounts(f_xy,f_x,f_y,N) )
}


#' Helper function to calculate probability from confusion matrix stats in dplyr friendly manner
#'
#' The purpose of this is to calculate the probabilities of binary outcomes from
#' a table of true pos, true neg, false pos and false neg trials 
#' 
#' @param df a dataframe containing one observation per row
#' @param tpVar the datatable column containing frequency of cooccurrence of true positives
#' @param fpVar the datatable column containing frequency of occurrence of false positives
#' @param fnVar the datatable column containing frequency of occurrence of false negatives
#' @param tnVar the datatable column containing frequency of occurrence of true negatives
#' @return the datatable with additional columns for all the probabilities associated with each outcome (i.e. a 2x2 confusion matrix)
#' @import dplyr
#' @export
probabilitiesFromConfusionMatrix = function(df, tpVar, fpVar, fnVar, tnVar) {
  tpVar = ensym(tpVar)
  fpVar = ensym(fpVar)
  fnVar = ensym(fnVar)
  tnVar = ensym(tnVar)
  return(
    df %>% mutate(
      total = !!tpVar+!!fpVar+!!fnVar+!!tnVar,
      p_x1y1 = as.double(!!tpVar)/total,
      p_x1y0 = as.double(!!fnVar)/total,
      p_x0y1 = as.double(!!fpVar)/total,
      p_x0y0 = as.double(!!tnVar)/total,
      p_x1 = p_x1y1+px1y0,
      p_x0 = p_x0y1+px0y0,
      p_y1 = p_x1y1+px0y1,
      p_y0 = p_x1y0+px0y0
    )
  )
}

#' Helper function to calculate probability from counts in dplyr friendly manner
#'
#' The purpose of this is to calculate the probabilities of given all binary outcomes from
#' a table of marginal frequencies in a dplyr friendly way. 
#' 
#' @param df a dataframe containing one observation per row
#' @param x1y1Var the datatable column containing frequency of cooccurrences of events x1 and y1
#' @param x1Var the datatable column containing frequency of occurrences of event x1
#' @param y1Var the datatable column containing frequency of occurrences of event y1
#' @param totalVar the datatable column containing total number of events
#' @return the datatable with additional columns for all the probabilities associated with each outcome (i.e. a 2x2 confusion matrix)
#' @import dplyr
#' @export
probabilitiesFromCounts = function(df, x1y1Var, x1Var, y1Var, totalVar) {
  x1y1Var = ensym(x1y1Var)
  x1Var = ensym(x1Var)
  y1Var = ensym(y1Var)
  totalVar = ensym(totalVar)
  return(
    df %>% mutate(
      p_x1 = as.double(!!x1Var)/!!totalVar,
      p_x0 = 1-p_x1,
      p_y1 = as.double(!!y1Var)/!!totalVar,
      p_y0 = 1-p_y1,
      p_x1y1 = as.double(!!x1y1Var)/!!totalVar,
      p_x1y0 = p_x1 - p_x1y1,
      p_x0y1 = p_y1 - p_x1y1,
      p_x0y0 = 1.0-p_x1y1-p_x0y1-p_x1y0
    ) %>% mutate(
      p_x1 = ifelse(p_x1 < 0.0 | p_x1 > 1.0, NA, p_x1),
      p_x0 = ifelse(p_x0 < 0.0 | p_x0 > 1.0, NA, p_x0),
      p_y1 = ifelse(p_y1 < 0.0 | p_y1 > 1.0, NA, p_y1),
      p_y0 = ifelse(p_y0 < 0.0 | p_y0 > 1.0, NA, p_y0),
      p_x1y1 = ifelse(p_x1y1 < 0.0 | p_x1y1 > 1.0, NA, p_x1y1),
      p_x0y1 = ifelse(p_x0y1 < 0.0 | p_x0y1 > 1.0, NA, p_x0y1),
      p_x1y0 = ifelse(p_x1y0 < 0.0 | p_x1y0 > 1.0, NA, p_x1y0),
      p_x0y0 = ifelse(p_x0y0 < 0.0 | p_x0y0 > 1.0, NA, p_x0y0)
    ) 
  )
}
