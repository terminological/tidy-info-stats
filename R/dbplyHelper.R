# N.B. test cases exist

#### DBPlyr grouped data utilities ----

#' Calculates a join list
#' 
#' @param df - a df which may be grouped
#' @param groupVars - the grouping for which we want to create a label as a list of columns quoted by vars(...)
#' @param defaultJoin - if there is no grouping we need one column to join by.
#' @export
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
    if (length(grps) == 0) stop("groupMutate can only work on grouped data")
    tmp_df <- df %>% group_by(!!!grps) %>% summarize(...)
    return(df %>% left_join(tmp_df, by=sapply(grps,as.character)))
  } else {
    return(df %>% mutate(...))
  }
}

#' Renames groups extending accross multiple columns do that subgroups have sequentail id
#' 
#' Sometime simpler to rename multiple grouping columns as one single column which encodes all the possible combinations.
#'  
#' @param df - a df which may be grouped in which case the grouping can be interpreted as a feature and each group must have the same number in it.
#' @param groupVars - the columns(s) for which we want to create a label - as a list of columns quoted by vars(...).
#' @param outputVar - the name of the target column for the new discrete variable (or sample identifier)
#' @param summarise - return dataframe as-is with additional column (FALSE - the default) or return dataframe as group summary with only grouping info and output (TRUE)
#' @param consistentBetweenGroups - if the groups represent discrete variables that vary between groups then they are not comparable. 
#' @return a dbplyr dataframe containing the grouped function
#' @export
labelGroup = function(df, groupVars, outputVar, summarise=FALSE, consistentBetweenGroups=FALSE) {
  outputVar = ensym(outputVar)
  grps = df %>% groups()

  if (consistentBetweenGroups) {
    
    tmp = df %>% select(!!!groupVars) %>% distinct() %>%
      arrange(!!!groupVars) %>% mutate(!!outputVar := row_number()) 
    
    if (summarise) {
      return(tmp)
    } else {
      joinList = df %>% ungroup() %>% joinList(groupVars)
      return(df %>% left_join(tmp,by=joinList))
    }
    
  } else {
    
    tmp = df %>% select(!!!grps,!!!groupVars) %>% distinct() %>%
      group_by(!!!grps) %>% arrange(!!!groupVars) %>% mutate(!!outputVar := row_number()) 
    
    if (summarise) {
      return(tmp)
    } else {
      joinList = df %>% joinList(groupVars)
      return(df %>% left_join(tmp,by=joinList))
    }
  }
}


#' Create an id based for each sample based on a unique combination of variables
#' 
#' Sometime simpler to rename multiple grouping columns as one single column which encodes all the possible combinations.
#'  
#' @param df - a df which may be grouped in which case the grouping can be interpreted as a feature and each group must have the same number in it.
#' @param groupVars - the columns(s) for which we want to create a label - as a list of columns quoted by vars(...). Essentially a combination of discrete variables. i.e. identifies each sample.
#' @param sampleIdVar - the name of the target column for the new discrete variable (or sample identifier)
#' @return a dbplyr dataframe containing the grouped function
#' @export
createSequentialIdentifier = function(df, groupVars, sampleIdVar) {
  sampleIdVar = ensym(sampleIdVar)
  grps = df %>% groups()
  
  tmp = df %>% select(!!!groupVars) %>% distinct() %>% arrange(!!!groupVars) %>% mutate(!!sampleIdVar := row_number()) 
  
  joinList = df %>% ungroup() %>% joinList(groupVars)
  return(df %>% left_join(tmp,by=joinList))
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

#### Digamma calculation ----

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


#### Savitsky Golay filter ----

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
  # TODO: some form of support range for each group (currently fixed from 0 to 1)
  # k_05=2
  # df = tibble(group = sample.int(4,size=100,replace=TRUE), value=rnorm(100)) %>% group_by(group)
  # continuousVar = "value"
  grps = df %>% groups()
  continuousVar = ensym(continuousVar)
  outputVar = ensym(outputVar)
  
  if ("tbl_sql" %in% class(df)) {
    
    # a dbplyr table - copy sgolay coefficients to database
    df = df %>% mutate(
      tmp_id=row_number(),
      N = n(),
      tmp_value = !!continuousVar,
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
        sample = tmp_id-selector) %>%
      mutate(
        position_match = ifelse(sample<=k_05, sample-k_05-1, ifelse(sample>(N-k_05),k_05-(N-sample),0))
      ) %>% filter(position == position_match) %>% compute()
    
    # coefficients are groupwise adjusted for sample size which may vary between group
    tmp = tmp %>% coefficientAdj(N, derivative = m)
    
    tmp2 = tmp %>% group_by(!!!grps,sample) %>% summarise(
      N = max(N,na.rm = TRUE), 
      k_05 = max(k_05,na.rm = TRUE),
      !!outputVar := sum(coefficientAdj*tmp_value, na.rm = TRUE), 
      !!continuousVar := sum(tmp_value*ifelse(tmp_id==sample,1.0,0.0), na.rm = TRUE)
    ) %>% mutate(
      !!outputVar := ifelse(N<(k_05*2+1),NA,!!outputVar)
    ) %>% compute()
    
    return(tmp2 %>% select(-sample, -k_05))
    
  } else {
    
    # Not a dplyr table - use signal::sgolayfilt
    out = df %>% rename(tmp_x_continuous = !!continuousVar) %>%  group_by(!!!grps) %>% group_modify(
      function(d,...) {
        k = k_05*2+1
        samples = nrow(d)
        # if (k >= samples) k = samples-samples%%2-1
        if (k < samples-1) {
          temp_est = signal::sgolayfilt(d$tmp_x_continuous, p=p, n=k, m=m, ts=1.0/samples)
          return(
            tibble(
              N = samples,
              !!continuousVar := d$tmp_x_continuous,
              !!outputVar := temp_est
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

#### Dataframe collection ----

#' Converts a tidy dataframe into a Matrix::sparseMatrix 
#' 
#' offloading the majority of processing onto sql if dbplyr tables are involved
#' 
#' @param df - a df
#' @param rowVar - the dataframe columns(s) which define the matrix row, quoted by vars(...) - typically this is the observation id
#' @param colVar - the dataframe columns(s) which define the matrix columns, quoted by vars(...) - typically this is the feature id
#' @param valueVar - the name of the value variable. (#TODO could be missing - in which case use binary)
#' @param rowNameVar - (optional) the dataframe column continaing the row names otherwise use rowVar
#' @param colNameVar - (optional) the dataframe column continaing the column names otherwise use colVar
#' @param ... - other parameters passes to Matrix::sparseMatrix
#' @return a dbplyr dataframe containing the grouped function
#' @export
collectAsSparseMatrix = function(df, rowVar, colVar, valueVar=NULL, rowNameVar=NULL, colNameVar=NULL, ...) {
  valueVar = tryCatch(ensym(valueVar), error=function(e) NULL)
  rowVar = ensym(rowVar)
  colVar = ensym(colVar)
  rowNameVar = tryCatch(ensym(rowNameVar), error=function(e) NULL)
  colNameVar = tryCatch(ensym(colNameVar), error=function(e) NULL)
  df = df %>% ungroup()
  
  # group data by rowVar and generate a sequential row_id and label
  if (identical(rowNameVar,NULL)) {
    rows = df %>% select(!!rowVar) %>% distinct() %>% arrange(!!rowVar) %>% mutate(tmp_row_id = row_number())
    rowLabels = rows %>% pull(!!rowVar)
  } else {
    rows = df %>% group_by(!!rowVar) %>% summarise(tmp_label = !!rowNameVar) %>% arrange(!!rowVar) %>% mutate(tmp_row_id = row_number())
    rowLabels = rows %>% pull(tmp_label) %>% make.unique(sep=" ")
  }
  
  # group data by colVar and generate a sequential col_id and label
  if (identical(colNameVar,NULL)) {
    cols = df %>% select(!!colVar) %>% distinct() %>% arrange(!!colVar) %>% mutate(tmp_col_id = row_number())
    colLabels = cols %>% pull(!!colVar)
  } else {
    cols = df %>% group_by(!!colVar) %>% summarise(tmp_label = min(!!colNameVar)) %>% arrange(!!colVar) %>% mutate(tmp_col_id = row_number())
    colLabels = cols %>% pull(tmp_label) %>% make.unique(sep=" ")
  }
  
  if (identical(valueVar,NULL)) {
    
    # There is no value so we use the pattern matrix form of Matrix::sparseMatrix
    data = df %>% 
      select(!!rowVar,!!colVar) %>% 
      inner_join(rows, by=as.character(rowVar)) %>% 
      inner_join(cols, by=as.character(colVar)) %>% 
      select(tmp_row_id, tmp_col_id) %>%
      collect()
    #browser()
    return(Matrix::sparseMatrix(i = data$tmp_row_id, j=data$tmp_col_id, dimnames= list(rowLabels,colLabels), ...))
    
  } else {
    
    data = df %>% 
      select(!!rowVar,!!colVar,!!valueVar) %>% 
      inner_join(rows, by=as.character(rowVar)) %>% 
      inner_join(cols, by=as.character(colVar)) %>% 
      rename(tmp_value = !!valueVar) %>%
      select(tmp_row_id, tmp_col_id, tmp_value) %>%
      collect()
    #browser()
    return(Matrix::sparseMatrix(i = data$tmp_row_id, j=data$tmp_col_id, x=data$tmp_value, dimnames=list(rowLabels,colLabels), ...))
    
  }
  
}



#' Converts a tidy dataframe into a Matrix::sparseMatrix 
#' 
#' offloading the majority of processing onto sql if dbplyr tables are involved
#' 
#' @param df - a df
#' @param rowVar - the dataframe columns(s) which define the matrix row, quoted by vars(...) - typically this is the observation id
#' @param outcomeVar - (optional) the dataframe column continaing the outcome
#' @param factorise - convert outcomeVar to a factor (default FALSE)
#' @param ... - other parameters passed to as.factor
#' @return a dbplyr dataframe containing the grouped function
#' @export
collectOutcomeVector = function(df, rowVar, outcomeVar, factorise = TRUE, ...) {
  rowVar = ensym(rowVar)
  outcomeVar = ensym(outcomeVar)
  df = df %>% ungroup()
  
  rows = df %>% select(!!rowVar, !!outcomeVar) %>% distinct() %>% arrange(!!rowVar)
  
  outcomeVector = rows %>% pull(!!outcomeVar)
  if (factorise) {
    outcomeVector = as.factor(outcomeVector, ...)
  }
  
  return(outcomeVector)
}




#' Converts a tidy dataframe into a Matrix::sparseMatrix of features & associated outcome vector
#' 
#' offloading the majority of processing onto sql if dbplyr tables are involved.
#' 
#' @param df - a df
#' @param sampleVar - the dataframe columns(s) which define the matrix row, quoted by vars(...) - typically this is the observation id ( see createSequentialIdentifier(...) )
#' @param outcomeVar - the dataframe columns(s) which define the matrix columns, quoted by vars(...) - typically this is the feature id
#' @param featureVar - the dataframe columns(s) which define the matrix columns, quoted by vars(...) - typically this is the feature id
#' @param valueVar - the name of the value variable. (#TODO could be missing - in which case use binary)
#' @param ... - other parameters passes to Matrix::sparseMatrix & as.factor (for outcomes)
#' @return a list with the following elements:
#' 
#' * rowLabels - the labels for each row in order
#' * colLabels - the feature labels in order
#' * matrix - a Matrix::sparseMatrix of the data, with values as doubles
#' * outcome - a vector of outcomes (probably as a factor)
#' @export
collectAsTrainingSet = function(df, sampleVar, outcomeVar, featureVar, valueVar=NULL, featureNameVar=NULL, factorise = TRUE, ...) {
  valueVar = tryCatch(ensym(valueVar), error=function(e) NULL)
  featureNameVar = tryCatch(ensym(featureNameVar), error=function(e) NULL)
  sampleVar = ensym(sampleVar)
  featureVar = ensym(featureVar)
  outcomeVar = ensym(outcomeVar)
  df = df %>% ungroup()
  
  out = list()
  
  # group data by rowVar and generate a sequential row_id and label
  rows = df %>% select(!!sampleVar, !!outcomeVar) %>% distinct() %>% arrange(!!sampleVar) %>% mutate(sampleIndex = row_number())
  
  nonUnique = rows %>% group_by(!!sampleVar) %>% summarise(count = n()) %>% filter(count > 1) %>% collect();
  if(nrow(nonUnique) != 0) {
    print(nonUnique)
    stop("non unique outcomes detected in data set")
  }
  
  out$rowLabels = rows %>% pull(!!sampleVar)
  out$outcome = rows %>% pull(!!outcomeVar)
  if (factorise) {out$outcome = as.factor(out$outcome, ...)}
  
  # group data by colVar and generate a sequential col_id and label
  if (identical(featureNameVar,NULL)) {
    cols = df %>% select(!!featureVar) %>% distinct() %>% arrange(!!featureVar) %>% mutate(featureIndex = row_number())
    out$colLabels = cols %>% pull(!!featureVar)
  } else {
    cols = df %>% group_by(!!featureVar) %>% summarise(featureLabel = min(!!featureNameVar)) %>% arrange(!!featureVar) %>% mutate(featureIndex = row_number())
    out$colLabels = cols %>% pull(featureLabel) %>% make.unique(sep=" ")
  }
  
  if (identical(valueVar,NULL)) {
    
    # There is no value so we use the pattern matrix form of Matrix::sparseMatrix
    data = df %>% 
      select(!!sampleVar,!!featureVar) %>% 
      inner_join(rows, by=as.character(sampleVar)) %>% 
      inner_join(cols, by=as.character(featureVar)) %>% 
      select(sampleIndex, featureIndex) %>%
      collect()
    #browser()
    out$matrix = Matrix::sparseMatrix(i = data$sampleIndex, j=data$featureIndex, dimnames= list(rowLabels,colLabels), ...)
    
  } else {
    
    data = df %>% 
      select(!!sampleVar,!!featureVar,!!valueVar) %>% 
      inner_join(rows, by=as.character(sampleVar)) %>% 
      inner_join(cols, by=as.character(featureVar)) %>% 
      rename(tmp_value = !!valueVar) %>%
      select(sampleIndex, featureIndex, tmp_value) %>%
      collect()
    #browser()
    out$matrix = Matrix::sparseMatrix(i = data$sampleIndex, j=data$featureIndex, x=data$tmp_value, dimnames=list(out$rowLabels,out$colLabels), ...)
    
  }

  return(out)
}

#' Sparse matrix conversion for liblinear
#' 
#' Converts a Matrix::sparseMatrix to a SparseM::matrix.csr without blowing up the computer
#' 
#' @param X as Matrix::sparseMatrix
#' @return a SparseM::matrix.csr
#' @export
sparseMatrixToSparseMCsr = function(X) {
  X.csc <- new("matrix.csc", ra = X@x,
               ja = X@i + 1L,
               ia = X@p + 1L,
               dimension = X@Dim)
  X.csr <- SparseM::as.matrix.csr(X.csc)
  return (X.csr)
}