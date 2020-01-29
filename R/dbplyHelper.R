#' Execute a mutate function on grouped data in all databases
#' 
#' Some databases don't support window functions over grouped data. This requires a workaround to group and summarise then join the data.
#' All database seem to support group_by(...) and then mutate(x = sum(...)) but not all will do group_by(...) and then mutate(x = mean(...)) for example.
#' 
#' Ranking functions all seem to work natively (including row_number)
#' 
#' * sum
#' * row_number
#' * rank
#' * lead, lag
#' 
#' Aggregation on the other hand functions are variably supported by the default mutate depending on the database
#' 
#' * mean
#' * sd
#' * min
#' * max
#' 
#' groupMutate allows any window functions to be applied to any database but using an intermediate table
#' 
#' @param df - a df which may be a dbplyr table
#' @param ... - the contents of the mutuate function
#' @return a dbplyr dataframe containing the grouped function
#' @export
groupMutate = function(df, ...) {
  if ("tbl_sql" %in% class(df)) {
    grps = df %>% groups()
    if (length(grps == 0)) stop("groupMutate can only work on grouped data")
    tmp_df <- df %>% group_by(!!!grps) %>% summarize(...)
    return(df %>% left_join(tmp_df, by=sapply(grps,as.character)))
  } else {
    return(df %>% mutate(...))
  }
}

#' Renames groups extending accross multiple columns
#' 
#' Sometime simpler to rename multiple grouping columns as one single column which encodes all the possible combinations.
#'  
#' @param df - a df which may be grouped
#' @param groupVars - the columns(s) for which we want to create a label - as a list of columns quoted by vars(...)
#' @param outputVar - the name of the target column for the new discrete variable
#' @param summarise - return dataframe as-is with additional column (FALSE - the default) or return dataframe as group summary with only grouping info and output (TRUE)
#' @return a dbplyr dataframe containing the grouped function
#' @export
labelGroup = function(df, groupVars, outputVar, summarise=FALSE) {
  outputVar = ensym(outputVar)
  grps = df %>% groups()

  tmp = df %>% select(!!!grps,!!!groupVars) %>% distinct() %>%
    group_by(!!!grps) %>% arrange(!!!groupVars) %>% mutate(!!outputVar := row_number()) 

  if (summarise) {
    return(tmp)
  } else {
    joinList = df %>% joinList(groupVars)
    return(df %>% left_join(tmp,by=joinList))
  }
  
}


#' Groupwise count
#' 
#' Performs a 2 level count either preserving the structure of the dataframe (e.g. as a mutate function) or as a summary, returning
#' the dataframe with total counts for the grouping (N) and for the subgroup defined by "groupVars" (N_x). This also allows us to 
#' calculate a probability of the subgroup in the group. This is different from an rank_percent in that the input datafram may have already been
#' summarised
#' 
#' @param df - a df which may be grouped. Grouping typically will be on a feature. N is the count of the items in the group
#' @param groupVars - the grouping for which we want to create a label as a list of columns quoted by vars(...). This could be an outcome and 
#' @param countVar - optional: the datatable column containing the observed frequency of the event X. If this is missing the row count will be used instead (i.e. assumes each row is an observation).
#' @param summarise - return dataframe as-is with additional columns (N, N_x, p_x) (FALSE - the default) or return dataframe as group summary with only grouping info and output (TRUE)
#' @return the grouped dataframe containing at a minumum, the df grouping columns, the groupVars columns, and a groupwise count of both levels of grouping labelled N and N_x and the groupwise p_x.
#' @export
#' @examples 
#' mtcars %>% group_by(cyl) %>% groupwiseCount(vars(gear), summarise=TRUE) %>% mutate(p_x = N_x/N)
groupwiseCount = function(df, groupVars, countVar=NULL, summarise=FALSE) {
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  grps = df %>% groups()
  
  if (summarise) {
    if (identical(countVar,NULL)) {
      # group count + mutate surm
      return(df %>% group_by(!!!grps, !!!groupVars) %>% summarise(N_x = n()) %>% ungroup() %>% group_by(!!!grps) %>% mutate(N = sum(N_x)))
    } else {
      # group sum + mutate surm
      return(df %>% group_by(!!!grps, !!!groupVars) %>% summarise(N_x = sum(!!countVar)) %>% ungroup() %>% group_by(!!!grps) %>% mutate(N = sum(N_x)))
    }
  } else {
    if (identical(countVar,NULL)) {
      # mutate count + mutate count
      return(df %>% group_by(!!!grps, !!!groupVars) %>% mutate(N_x = n()) %>% ungroup() %>% group_by(!!!grps) %>% mutate(N = n()))
    } else {
      # group count + mutate sum + join to original
      joinList = df %>% joinList(groupVars)
      tmp = df %>% group_by(!!!grps, !!!groupVars) %>% summarise(N_x = sum(!!countVar)) %>% ungroup() %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
      return(df %>% left_join(tmp,by=joinList))
    }
  }
}


#' Check table is from a database and enforce collection
#' 
#' Fail if table is a dbplyr table and the user does not request collection
#' 
#' @param df - a dataframe
#' @param collect - boolean, should the table be collected if it is a dbplyr table.
#' @return the dataframe, collected locally.
#' @export
collectDf = function(df, collect) {
  if ("tbl_sql" %in% class(df)) {
    if (collect) {
      return(df %>% collect())
    } else {
      stop("This implementation does not support dbplyr data frames. use option collect=TRUE to coerce to a local data frame.")
    }
  }
  return(df)
}

#' a set of the first 50 digamma values
digammaLookup = tibble(n = c(1:50), digamma = digamma(c(1:50)))

#' SQL Digamma function
#' 
#' Calculate an estimate of the digamma function value of integers for dplyr and dbplyr dataframes. This will be estimated in SQL if needed.
#' 
#' @param df - a dataframe
#' @param inputVar - the reference to the column for which the digamma value will be calculated (this must be an integer)
#' @param outputVar - the column that the digamma result will be written to
#' @return the dataframe with an "outputVar" column containing the digamma value
#' @export
calculateDigamma = function(df, inputVar, outputVar) {
  inputVar = ensym(inputVar)
  outputVar = ensym(outputVar)
  if ("tbl_sql" %in% class(df)) {
    tmp = df %>% 
      left_join(digammaLookup %>% rename(!!inputVar := n, !!outputVar := digamma), by=as.character(inputVar), copy=TRUE) %>%
      mutate(!!outputVar := ifelse(is.na(!!outputVar), log(!!inputVar-1) + 1.0/(2*(!!inputVar-1) - 1.0/(12*(!!inputVar-1)^2) ), !!outputVar))
    return(tmp)
  } else {
    return(df %>% mutate(!!outputVar := digamma(!!inputVar)))
  }
}

#' Calculates a join list
#' 
#' @param df - a df which may be grouped
#' @param groupVars - the grouping for which we want to create a label as a list of columns quoted by vars(...)
#' @param defaultJoin - if there is no grouping we need one column to join by.
joinList = function(df,groupVars=NULL,defaultJoin=NULL) {
  grps = df %>% groups()
  joinList = c()
  if (!identical(defaultJoin,NULL)) {
    joinList = c(joinList,defaultJoin)
  }
  if (length(grps)!=0) {
    joinList = c(joinList, sapply(grps,as.character))
  }
  if (!identical(groupVars,NULL)) {
    joinList = c(joinList, as.vector(sapply(groupVars,as_label)))
  }
  return(joinList)
}

#' Precalculates a sgolay filter (WIP)
#'
#' @param polynomialOrder - order of the poynomial to fit
#' @param filterLength - length of the filter
#' @param derivative - derivative order for the fitted function
#' @return a tidy format sgolay coefficient dataframe where position
#' @export
#' @examples 
#' sgolayTable(2,5,1)
sgolayTable = function(polynomialOrder, filterLength, derivative) {
  return(tibble(
    coefficient = as.vector(signal::sgolay(p=polynomialOrder, n=filterLength, m=derivative ,1)),
    position = rep(seq(-floor(filterLength/2),floor(filterLength/2)),filterLength),
    offset = as.vector(sapply(seq(-floor(filterLength/2),floor(filterLength/2)),rep,filterLength))
  ))
}

#' adjust the baseline coefficients of a SGoly filter for specifc derivative order and sample size (WIP)
#' @param sgolayTable - a data frame containing the sgolay coefficients joined to a data set
#' @param sampleSizeVar - the colum containing the sample size it is to be applied to
#' @param supportRange - the total support range of the data being smoothed (defaults to 1 for CDFs)
#' @export
coefficientAdj = function(sgolayTable, sampleSizeVar, derivative, supportRange=1) {
  sampleSizeVar = ensym(sampleSizeVar)
  sgolayTable %>% mutate(coefficientAdj = coefficient / ((as.double(supportRange)/(as.double(!!sampleSizeVar)))^local(derivative)) )
}

# tibble(group = c(1,1,1,2,2,2,2), value=runif(7)) %>% group_by(group) %>% mutate(N=n())
# tibble(group = sample.int(4,size=100,replace=TRUE), value=rnorm(100))

#' apply a sgolay filter groupwise to data in a dbplyr dataframe in SQL friendly manner
#' @param df - a potentially grouped data frame with continuous observations
#' @param continuousVar - the column containing the samples
#' @param outputVar - the column to write the filter result
#' @param k_05 - the half value of the filter width
#' @param p - the filter order
#' @param m - the mth derivative
#' @return a summarised dataframe including group info, the continuous variable and an output variable representing the filtered value
#' @export
applySGolayFilter = function(df, continuousVar, outputVar, k_05, p, m) {
  # k_05=2
  # df = tibble(group = sample.int(4,size=100,replace=TRUE), value=rnorm(100)) %>% group_by(group)
  # continuousVar = "value"
  grps = df %>% groups()
  continuousVar = ensym(continuousVar)
  outputVar = ensym(outputVar)
  
  if ("tbl_sql" %in% class(df)) {
    
    # a dbplyr table - copy sgolay coefficients to database
    df = df %>% arrange(!!continuousVar) %>% mutate(
      tmp_id=row_number(),
      N = n(),
      k_05 = k_05 # ifelse(N<7,NA,ifelse(k_05<floor((N-1)/2),k_05,floor((N-1)/2)))
    ) %>% mutate(
      maxSelector = ifelse(tmp_id<=(k_05*2+1),tmp_id-1,k_05),
      minSelector = ifelse(tmp_id>=(N-k_05*2),tmp_id-N,-k_05),
      # offset = ifelse(tmp_id<=k_05, tmp_id-k_05-1, ifelse(tmp_id>(N-k_05),k_05-(N-tmp_id),0)),
      tmp_join=1
    ) 
    
    # the coefficients are calculated assuming a ts of 1. This is later adjusted for sample size
    coeff = sgolayTable(p, k_05*2+1, m) %>% mutate(tmp_join=1)
    coeff = coeff %>% mutate(selector = offset-position)
    
    tmp = df %>% left_join(coeff, by="tmp_join", copy=TRUE) %>% select(-tmp_join) %>% 
      filter(maxSelector>=selector & minSelector<=selector) %>% 
      mutate(
        sample = tmp_id-selector,
        position_match = ifelse(sample<=k_05, sample-k_05-1, ifelse(sample>(N-k_05),k_05-(N-sample),0))
      ) %>% filter(position == position_match) 
    
    # coefficients are groupwise adjusted for sample size which may vary between group
    tmp = tmp %>% coefficientAdj(N, derivative = m)
    
    
    tmp2 = tmp %>% group_by(!!!grps,sample) %>% summarise(
      N = max(N,na.rm = TRUE), 
      k_05 = max(k_05,na.rm = TRUE),
      !!outputVar := sum(coefficientAdj*value), 
      !!continuousVar := sum(value*ifelse(tmp_id==sample,1.0,0.0))
    ) %>% mutate(
      !!outputVar := ifelse(N<(k_05*2+1),NA,!!outputVar)
    )
    
    return(tmp2 %>% select(-sample, -k_05))
    
  } else {
    
    # Not a dplyr table - use signal::sgolayfilt
    out = df %>% rename(tmp_x_continuous = !!continuousVar) %>%  group_by(!!!grps) %>% arrange(tmp_x_continuous) %>% group_modify(
      function(d,...) {
        k = k_05*2+1
        samples = nrow(d)
        # if (k >= samples) k = samples-samples%%2-1
        if (k < samples-1) {
          d_x_d_r = signal::sgolayfilt(d$tmp_x_continuous, p=p, n=k, m=m, ts=1.0/samples)
          return(
            tibble(
              N = samples,
              !!continuousVar := d$tmp_x_continuous,
              !!outputVar := d_x_d_r # prevent negative gradient - d_x_d_r is an inverse cdf - always positive
            )
          )
        } else {
          return(
            tibble(
              N = samples, # TODO
              !!continuousVar := d$tmp_x_continuous,
              !!outputVar := rep(NA,length(d$tmp_x_continuous))
            )
          )
        }
      })
      return(out)
  }
  
  
  
}