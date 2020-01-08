#' calculate entropy of a sequence of 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param continuousVar - the columns that define the discrete subgroups of the data.
#' @param discreteOutputVar - the name of the value to create in the dataframe
#' @param method - the method employed - valid options are "ByRank", "ByValue", "Manual"
#' @param ... - the other parameters are passed onto the implementations
#' @return a single value for the entropy of the vector
#' @import dplyr
#' @export
discretise = function(df, continuousVar, discreteOutputVar, method, ...) {
  switch (method,
          ByRank = discretise_ByRank(df, {{continuousVar}}, {{discreteOutputVar}}, ...),
          ByValue = discretise_ByValue(df, {{continuousVar}}, {{discreteOutputVar}}, ...),
          Manual = discretise_Manual(df, {{continuousVar}}, {{discreteOutputVar}}, ...),
          #Histogram = discretise_Histogram(df, {{continuousVar}}, {{discreteOutputVar}}, ...)
          #ChiSq
          #MaxMI
          {stop(paste0(method," not a valid option"))}
  )
}

#TODO: non linear binning strategies

#' A binning strategy generator. Creates a function that can be used to create discretisation cuts based on the statistical paramaters of the data to be discretised
#' in particular the mean and sd of the observed data (in a group by group fashion)
#' 
#' @param bins - the number of bins to create out of the centiles of a log normal distribution
logNormalCentiles = function(bins) {
  return(
    function(mean,sd,...) {
      # TODO:
      # need to solve these two equations for μ & σ
      # mean = exp(μ + 1/2 σ^2)
      # sd = sqrt(exp(2*μ + σ^2)*(exp(σ^2) - 1))
      # use DistributionFunctions
      return(tibble(
        cut = qlnorm(seq(0,1,length.out = bins+1), meanlog = log(mean), sdlog = sd)[2:bins]
      ))
    }
  )
}

#' A binning strategy generator. Creates a function that can be used to create discretisation cuts based on the statistical paramaters of the data to be discretised
#' in particular the minimum and maximum value of the observed data (in a group by group fashion)
#' 
#' @param bins - the number of bins to create for each group
fixedNumber = function(bins) {
  return(
    function(n,min,max, ...) {
      return(tibble(
        cut = seq(min,max,length.out = bins+1)[2:bins]
      ))
    }
  )
}

#' A binning strategy generator. Creates a function that can be used to create discretisation cuts based on the statistical paramaters of the data to be discretised
#' in particular the sample size, minimum and maximum value of the observed data (in a group by group fashion)
#' 
#' @param slope - the number of observations for each bin if the observation numbers fall between minBins*slope and maxBins*slope
#' @param minBins - the smallest number of bins to create for each group
#' @param maxBins - the largest number of bins to create for each group
#' @param fn - the distribution function that genertes the cut points (e.g. tidyinfostats::fixedNumber (uniform), tidyinfostats::logNormalCentiles  )
linearBySize = function(slope,minBins,maxBins, fn=fixedNumber) {
  return(
    function(n, ...) {
      # this function will return a tibble of cut points given the n in group, min in group, max in group
      bins = as.integer(ifelse(n<minBins*slope,minBins,ifelse(n>maxBins*slope,maxBins,as.double(n)/slope)))
      return(fn(bins)(n, ...))
    }
  )
}

#' A binning strategy generator. Creates a function that can be used to create discretisation cuts based on the statistical paramaters of the data to be discretised
#' in particular the minimum and maximum value of the observed data (in a group by group fashion)
#' 
#' @param width - the width of each bin
fixedWidth = function(width) {
  return(
    function(min, max, ...) {
      return(tibble(
        cut = seq(min,max,by=width)
      ))
    }
  )
}


#' Discretise data by the value of a continuous variable
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param continuousVar - the columns that define the continuous data.
#' @param discreteOutputVar - the name of the value to create in the dataframe for the discrete data
#' @param bins - (optional) number of bins
#' @param binStrategy - if the number of bins is not set they will be calculated using this bin strategy
#' @return a dataframe with an additional column with the discrete categories
#' @import dbplyr
#' @export
discretise_ByValue = function(df, continuousVar, discreteOutputVar, bins=NA, binStrategy=linearBySize(slope=8,minBins=4,maxBins=100), ...) {
  
  discreteOutputVar = ensym(discreteOutputVar)
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  if (!is.na(bins)) binStrategy = fixedNumber(bins)
  
  cutsDf = df %>% summarise(
    N=n(),
    min = min(!!continuousVar),
    max = max(!!continuousVar),
    mean = mean(!!continuousVar),
    sd = sd(!!continuousVar)
  ) %>% collect()
  
  cutsDf = cutsDf %>% group_by(!!!grps) %>% group_modify(function(d,...) { binStrategy(n=d$N,min=d$min,max=d$max,mean=d$mean,sd=d$sd) })
  
  tmp2 = discretise_Manual(df, !!continuousVar, !!discreteOutputVar, cutsDf, ...)
  
  return(tmp2)
  
}

#' Discretise data by the rank of a continuous variable
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param continuousVar - the columns that define the continuous data.
#' @param discreteOutputVar - the name of the value to create in the dataframe for the discrete data
#' @param bins - (optional) number of bins
#' @param binStrategy - if the number of bins is not set they will be calculated using this bin strategy
#' @param ... - other options passed onto tidyinfostats::discretise_Manual
#' @return a dataframe with an additional column with the discrete categories
#' @import dbplyr
#' @export
discretise_ByRank = function(df, continuousVar, discreteOutputVar, bins=NA, binStrategy=linearBySize(slope=8,min=4,max=100), ...) {
  
  # TODO: support factorise=TRUE - generate labels etc.
  
  discreteOutputVar = ensym(discreteOutputVar)
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  
  if (!is.na(bins)) {
    # the easy case - a fixed number of bins for all groups
    tmp = df %>% mutate(!!discreteOutputVar := ntile(!!continuousVar, n=bins) )
    # TODO: so if we want a label this is a bit harder
    # nead to get a group & tile wise max value and lag it for the other end
    return(tmp)
  } else {
    # we need to calculate a group wise binning cut-of based on rank
    tmp = df %>% mutate(tmp_rank = min_rank(!!continuousVar))
    
    cutsDf = tmp %>% summarise(
      N=n()
    ) %>% collect()
    
    cutsDf = cutsDf %>% group_by(!!!grps) %>% group_modify(function(d,...) { binStrategy(n=d$N,min=1,max=d$N) }) %>% mutate(cut=as.integer(cut))
    # unfortunately ntile cannot cut uneven groups in dbplyr setting so we do manually
    tmp2 = discretise_Manual(tmp, tmp_rank, !!discreteOutputVar, cutsDf, factorise=TRUE, ...)
    tmp2 = tmp2 %>% select(-tmp_rank)
    # convert meaningless label (based on rank) into integer factor - if we want a meaningful labelled factor we will have to create it from the result.
    # a max(continuousVar) from a group and discreteOutputVar wise grouping with a max value will get one end of the label
    tmp2 = tmp2 %>% mutate(!!discreteOutputVar := as.integer(!!discreteOutputVar))
    return(tmp2)
  }
  
}

#' generate a data fram of cuts for every group in a dataframe from a vector of cuts, and the dataframe
#' 
#' @param df - a grouped dataframe
#' @param cuts - a vector of cuts - not including the -Inf, Inf ends
#' @return a datafram of cuts
#' @export
cutsDfFromVector = function(df, cuts) {
  grps = df %>% groups()
  return(
    df %>% select(!!!grps) %>% distinct() %>% mutate(tmp_join=1) %>% left_join(tibble(tmp_join=1,cut=cuts),by="tmp_join") %>% select(-tmp_join)
  )
}

#' Discretise data using pre-defined break points
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param continuousVar - the columns that define the discrete subgroups of the data.
#' @param discreteOutputVar - the name of the value to create in the dataframe
#' @param cutsDf - manual break points as a dataframe containing cut points for every group (see cutsDfFromVector). 
#' @param lowerBounded - default FALSE, should the lower or upper values be included in the group
#' @param factorise - convert discrete values into an ordered factor (alternative is a character string). This is only useful if you have a single set of cuts for all groups.
#' @param noUnicode - by default unicode characters are not used for the label if the target is a dbplyr table
#' @return a dataframe containing the discreteOutputVar column
#' @export
discretise_Manual = function(df, continuousVar, discreteOutputVar, cutsDf, lowerBounded = FALSE, factorise = FALSE, noUnicode = ("tbl_sql" %in% class(df)), ...) {
  discreteOutputVar = ensym(discreteOutputVar)
  continuousVar = ensym(continuousVar)
  
  if (as.character(discreteOutputVar) %in% names(df)) stop(paste0("discrete variable already exists in df: ",as.character(discreteOutputVar)))
  
  grps = df %>% groups()
  if (length(grps)==0) {
    joinList = c("tmp_join")
  } else {
    joinList = c(sapply(grps,as.character), "tmp_join")
  }
  
  # cutsDf is a local R construct. it does not need to be dbplyr friendly
  # generate the ends of the data
  # TODO: split this out into a function and allow custom labelling
  cutsDf = cutsDf %>% select(!!!grps,cut) %>% group_by(!!!grps) %>% arrange(cut) %>% mutate(tmp_upper = cut, tmp_lower = lag(cut,1,NA)) %>% 
    union(cutsDf %>% group_by(!!!grps) %>% arrange(cut) %>% mutate(tmp_upper = lead(cut,1,NA), tmp_lower = cut)) %>% 
    mutate(tmp_join=1) %>% select(-cut) %>% distinct()
  
  if (noUnicode) {
    # cuts is an R dataframe but it is going to be used in an sql database
    # no unicode label for database
    cutsDf = cutsDf %>% mutate(!!discreteOutputVar := 
         ifelse(is.na(tmp_lower),paste0(ifelse(lowerBounded,"<","<="),tmp_upper),
                ifelse(is.na(tmp_upper),paste0(ifelse(lowerBounded,">=",">"),tmp_lower),
                       paste0(tmp_lower,"-",tmp_upper)
                )))
  } else {
    # label is unicode
    cutsDf = cutsDf %>% mutate(!!discreteOutputVar := 
           ifelse(is.na(tmp_lower),paste0(ifelse(lowerBounded,"<","\u2264"),tmp_upper),
                ifelse(is.na(tmp_upper),paste0(ifelse(lowerBounded,"\u2265",">"),tmp_lower),
                       paste0(tmp_lower,"\u2013",tmp_upper)
                )))
  }
  
  if (factorise) { 
    
    l = cutsDf %>% ungroup() %>% select(tmp_upper,!!discreteOutputVar) %>% distinct() %>% arrange(tmp_upper) %>% pull(!!discreteOutputVar)
    cutsDf = cutsDf %>% mutate(!!discreteOutputVar := factor(!!discreteOutputVar, ordered=TRUE, levels = l))
    
    if ("tbl_sql" %in% class(df)) {
      # discrete value converted to an integer for factor in database
      cutsDf = cutsDf %>% mutate(!!discreteOutputVar := as.integer(!!discreteOutputVar))
    }
    
  }
    
  # convert the data
  
  tmp = df %>% mutate(tmp_value = !!continuousVar) %>% mutate(tmp_join=1) # %>% select(-!!discreteOutputVar)
  if (lowerBounded) {
    tmp2 = tmp %>% left_join(cutsDf, by=joinList, copy=TRUE) %>% filter(is.na(tmp_lower) | tmp_lower<=tmp_value) %>% filter(is.na(tmp_upper) | tmp_upper>tmp_value) 
  } else {
    tmp2 = tmp %>% left_join(cutsDf, by=joinList, copy=TRUE) %>% filter(is.na(tmp_lower) | tmp_lower<tmp_value) %>% filter(is.na(tmp_upper) | tmp_upper>=tmp_value) 
  }
  tmp2 = tmp2 %>% select(-c(tmp_lower,tmp_upper,tmp_value,tmp_join))
  
  return(tmp2)
}

