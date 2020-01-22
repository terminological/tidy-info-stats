#' Helper function to calculate probability from grouped data in a tidy friendly manner
#'
#' The purpose of this is to calculate the probabilities of all binary outcomes from class data. 
#' This function is useful when you have 2 types of events (X and Y) and
#' you either have a set of observations of their co-occurrence, containing non-unique X & Y combinations, or you have a confusion matrix of the counts of their combinations where 
#' each row has a unique combination of X and Y and a third column contains the counts of XY co-occurrences.
#' 
#' @param df a dataframe containing 2 columns defining class of event X and class of event Y and either one row per event, 
#' or a count of observations, for each class combination. 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param groupXVars the datatable column(s) defining the class of event X quoted by vars(...) e.g. outcomes
#' @param groupYVars the datatable column(s) defining the class of event Y quoted by vars(...) e.g. predictions
#' @param countVar the datatable column containing the observed frequency combination of event XY. If this is missing the row count will be used instead
#' @return A new datatable with all possible combinations of X&Y and the probabilities associated with each outcome (i.e. an N(X) by N(Y) set of binary confusion matrices)
#' @import dplyr
#' @export
probabilitiesFromCooccurrence = function(df, groupXVars, groupYVars, countVar=NULL) {
  grps = df %>% groups()
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList,"tmp_join")
  if (identical(countVar,NULL)) {
    # there is no count column.  The data is grouped and has event X and event Y entried for each occurrence
    df = df %>% group_by(!!!grps, !!!groupXVars, !!!groupYVars) %>% summarise(N_xy = n())
  } else {
    # there is a count column.  The data is grouped and has event X and event Y and counts of occurrences
    df = df %>% group_by(!!!grps, !!!groupXVars, !!!groupYVars) %>% summarise(N_xy = sum(!!countVar))
  }
  N = df %>% ungroup() %>% group_by(!!!grps) %>% summarise(
    N=sum(N_xy)
    #, this was for miller-madow adj we can do this differently now in tidyinfostats
    # classes_XY=n_distinct(!!groupXVar,!!groupYVar),
    # classes_X=n_distinct(!!groupXVar), 
    # classes_Y=n_distinct(!!groupYVar)
    ) %>% 
    # mutate(mm_adjust = as.double(classes_X+classes_Y-classes_XY-1)/(2*N)) %>% 
    mutate(tmp_join=1)
  # N = N %>% mutate(classes_XY = nrow(N))
  X = df %>% ungroup() %>% group_by(!!!grps,!!!groupXVars) %>% summarise(N_x = sum(N_xy)) %>% mutate(tmp_join=1) # grouping
  # X = X %>% mutate(classes_X = nrow(X))
  Y = df %>% ungroup() %>% group_by(!!!grps,!!!groupYVars) %>% summarise(N_y = sum(N_xy)) %>% mutate(tmp_join=1) # grouping
  #Y = Y %>% mutate(classes_Y = nrow(Y))
  # I__mm = I__emp + (|X| + |Y| - |XY| - 1)/2N
  
  
  XY = (X %>% inner_join(Y, by=joinList) %>% inner_join(N, by=joinList)) %>% select(-tmp_join)
  joinAll = c(grpsList,as.vector(sapply(groupXVars,as_label)),as.vector(sapply(groupYVars,as_label)))
  XY = XY %>% 
    left_join(df, by=joinAll) %>% 
    mutate( N_xy = ifelse(is.na(N_xy),0,N_xy))
  return( XY %>% group_by(!!!grps) %>% probabilitiesFromCounts(N_xy,N_x,N_y,N) )
}

#' Helper function to calculate probability from discrete data in a tidy friendly manner
#'
#' The purpose of this is to calculate the probabilities of events from discrete data. 
#' This function is useful when you have either a set of observations of its occurrence, 
#' containing non-unique x events, or you have a counts of their events where 
#' each row has the type of observation of X=x and the countVar column contains the counts of the event.
#' 
#' @param df a dataframe containing columns defining class of observations of discrete variable X and either one row per observation, 
#' or a count of observations for each of the observed values of X. 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param groupVars the datatable column(s) defining the class of the observation quoted by vars(...)
#' @param countVar the datatable column containing the observed frequency of the event X. If this is missing the row count will be used instead.
#' @return A summary datatable with possible values of X and the total (N), the total count of that group (N_x) the probability (p_x), and self information (I_x) associated with the value of X
#' @import dplyr
#' @export
probabilitiesFromDiscrete = function(df, groupVars, countVar=NULL) {
  grps = df %>% groups()
  
  # groupwise count creates an N and N_x  column based on groupVars, and countVar
  df = df %>% groupwiseCount(groupVars, countVar, summarise=TRUE) %>% mutate(
    p_x = N_x/N,
    I_x = -log(p_x) #I_x is self information
    # entropy of X as an empirical measure will be sum(p_x*I_x) - the average of self information
  )
}



#' Helper function to calculate observed class counts of discrete variable X in grouped data in a tidy friendly manner. 
#' Class counts are needed for certain corrections to the entropy and 
#'
#' @param df a dataframe containing 2 columns defining class of event X and class of event Y and either one row per event, 
#' or a count of observations, for each class combination. 
#' df may also be grouped and in which case the grouping is preserved in the result.
#' @param groupVars the datatable column(s) defining the class of the observation quoted by vars(...)
#' @param countVar the datatable column containing the observed frequency combination of event XY. If this is missing the row count will be used instead
#' @param summarise - return a mutated (FALSE) or summarised (TRUE) result of the df
#' @return A mutated datatable with the count of possible values of X as C (repeated for every observation of X) and the corresponding maximum entropy of the source (max_H)
#' @import dplyr
#' @export
classCountFromGroup = function(df, groupVars, summarise=FALSE) {
  grps = df %>% groups()
  tmp = df %>% ungroup() %>% group_by(!!!grps) %>% groupMutate(
    C = n_distinct(!!!groupVars),
    max_H = log(C)
  ) 
  if (summarise) {
    return(tmp)
  } else {
    joinList = df %>% joinList(groupVars)
    return(df %>% left_join(tmp,by=joinList))
  }
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
      tmp_total = !!tpVar+!!fpVar+!!fnVar+!!tnVar,
      p_x1y1 = as.double(!!tpVar)/tmp_total,
      p_x1y0 = as.double(!!fnVar)/tmp_total,
      p_x0y1 = as.double(!!fpVar)/tmp_total,
      p_x0y0 = as.double(!!tnVar)/tmp_total,
      p_x1 = p_x1y1+px1y0,
      p_x0 = p_x0y1+px0y0,
      p_y1 = p_x1y1+px0y1,
      p_y0 = p_x1y0+px0y0
    ) %>% select(-tmp_total)
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
