#' calculate mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "KWindow","KNN","SGolay","DiscreteByRank","DiscreteByValue","Compression","Entropy"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI = function(df, discreteVars, continuousVar, method="KWindow", ...) {
  switch (method,
    KWindow = calculateDiscreteContinuousMI_KWindow(df, discreteVars, {{continuousVar}}, ...),
    KNN = calculateDiscreteContinuousMI_KNN(df, discreteVars, {{continuousVar}}, ...),
    SGolay = calculateDiscreteContinuousMI_SGolay(df, discreteVars, {{continuousVar}}, ...),
    Kernel = calculateDiscreteContinuousMI_Kernel(df, discreteVars, {{continuousVar}}, ...),
    DiscretiseByRank = calculateDiscreteContinuousMI_DiscretiseByRank(df, discreteVars, {{continuousVar}}, ...),
    DiscretiseByValue = calculateDiscreteContinuousMI_DiscretiseByValue(df, discreteVars, {{continuousVar}}, ...),
    Compression = calculateDiscreteContinuousMI_Compression(df, discreteVars, {{continuousVar}}, ...),
    Entropy = calculateDiscreteContinuousMI_Entropy(df, discreteVars, {{continuousVar}}, ...),
    {stop(paste0(method," not a valid option"))}
    
  )
}

#' summary info about the observations
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param expected - the total number of samples for each outcome expected if there were no missing values.
summariseObservations = function(df, discreteVars, countVar = NULL) {
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  grps = df %>% groups()
  outerJoin = df %>% joinList(discreteVars)
  
  left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
  right = df %>% ungroup() %>% select(!!!discreteVars) %>% distinct() %>% mutate(join=1)
  cross = left %>% left_join(right, by="join")
  
  cross = cross %>% left_join(df %>% groupwiseCount(discreteVars, !!countVar, summarise=TRUE), by=outerJoin) %>%
    mutate(
      N_x = ifelse(is.na(N_x), 0, N_x),
      N = ifelse(is.na(N), 0, N)
    )
  return(cross)
}

#' summary info
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param expected - the total number of samples for each outcome expected if there were no missing values.
expectedObservations = function(df, discreteVars, expected) {
  grps = df %>% groups()
  left = df %>% select(!!!grps) %>% distinct() %>% mutate(join=1)
  right = df %>% ungroup() %>% select(!!!discreteVars) %>% distinct() %>% mutate(join=1)
  cross = left %>% left_join(right, by="join")
  cross = cross %>% mutate(N_x = expected) %>% group_by(!!!grps) %>% mutate(N = sum(N_x))
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) where a total number of possible outcomes is known.
#' 
#' This corrects normal mutual information calculations for information carried by the absense of a variable. This is relevant for 
#' sparse data sets with many features such as NLP terms.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param expectedCount - a dataframe containing colums for grouping, and discreteVars, N and N_x columns with counts of
#' @param method - the method employed - valid options are "KWindow","KNN","SGolay","DiscreteByRank","DiscreteByValue","Compression"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousUnlabelledMI = function(df, discreteVars, continuousVar, expectedCount, method="KWindow", ...) {
  grps = df %>% groups()
  I_given_present = calculateDiscreteContinuousConditionalMI(df, discretevars, {{continuousVar}}, method, ...)
  outerJoinCols = df %>% joinList(discreteVars)
  tmp = expectedCount %>% left_join(I_given_present %>% rename(N_lab = N, N_x_lab = N_x, I_given_lab_and_x = I_x), by=outerJoinCols) %>% 
    mutate(
      N_unlab = N - N_lab,
      N_x_unlab = N_x - N_x_lab
    ) %>% mutate(
      # we need to know I_given_unlabelled_and_x
      p_lab = as.double(N_lab)/N,
      p_unlab = as.double(N_unlab)/N,
      p_unlab_given_x = as.double(N_x_unlab)/N_x,
      p_lab_given_x = as.double(N_x_lab)/N_x,
      I_given_unlab_and_x = ifelse(p_unlab_given_x ==0,0,-log(p_unlab_given_x))
    ) %>% mutate(
      # combine labelled and unlabelled I (independent)
      p_x = as.double(N_x)/N,
      I_given_x = p_lab_given_x*I_given_lab_and_x + p_unlab_given_x * I_given_unlab_and_x
    ) 
    
  tmp2 = tmp %>% group_by(!!!grps, method, N) %>% summarise(
    I = sum(p_x * I_given_x),
    I_sd = NA
  ) %>% mutate(
    method = paste0("Unlabelled - ",method)
  )
  
  return(tmp2)
}


#' calculate a pointwise mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "KWindow","KNN","Kernel","Quantile"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousConditionalMI = function(df, discreteVars, continuousVar, method="KWindow", ...) {
  switch (method,
          KWindow = calculateDiscreteContinuousConditionalMI_KWindow(df, discreteVars, {{continuousVar}}, ...),
          KNN = calculateDiscreteContinuousConditionalMI_KNN(df, discreteVars, {{continuousVar}}, ...),
          Kernel = calculateDiscreteContinuousConditionalMI_Kernel(df, discreteVars, {{continuousVar}}, ...),
          Entropy = calculateDiscreteContinuousConditionalMI_Entropy(df, discreteVars, {{continuousVar}}, ...),
          # DiscretiseByRank = calculateDiscreteContinuousMI_DiscretiseByRank(df, discreteVars, {{continuousVar}}, ...),
          # DiscretiseByValue = calculateDiscreteContinuousMI_DiscretiseByValue(df, discreteVars, {{continuousVar}}, ...),
          # Compression = calculateDiscreteContinuousMI_Compression(df, discreteVars, {{continuousVar}}, ...),
          {stop(paste0(method," not a valid option"))}
          
  )
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a kernel density estimator of PDF
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_Kernel = function(df, discreteVars, continuousVar, collect=FALSE, ...) { 
  grps = df %>% groups()
  
  tmp5 = calculateDiscreteContinuousConditionalMI_Kernel(df, discreteVars, {{continuousVar}}, collect, ...)
  
  tmp6 = tmp5 %>% group_by(!!!grps) %>% summarise(
    N = max(N),
    I = sum(p_x * I_given_x),
    I_sd = sum(p_x * I_given_x_sd),
    method = "Kernel"
  )
  
  return(tmp6)
}

#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y). I.e. the self information of X conditioned on each of the possible values for X
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousConditionalMI_Kernel = function(df, discreteVars, continuousVar, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  continuousVar = ensym(continuousVar)
  
  tmp = df %>% 
    mutate(y_continuous=!!continuousVar) %>%
    group_by(!!!grps) %>%
    mutate(y_min = min(y_continuous),y_max = max(y_continuous))
  
  tmp2 = tmp %>%  
    group_by(!!!grps,!!!discreteVars) %>% 
    probabilitiesFromContinuous_Kernel(y_continuous,y_min,y_max) %>%
    rename(N_x = N, p_y_given_x=p_x, I_y_given_x = I_x)
  
  tmp3 = tmp %>%
    group_by(!!!grps) %>%
    probabilitiesFromContinuous_Kernel(y_continuous,y_min,y_max) %>%
    rename(N = N, p_y = p_x, I_y = I_x)
  
  joinList = tmp %>% joinList(defaultJoin = "y_continuous")
    
  tmp4 = tmp2 %>% 
    inner_join(tmp3, by=joinList) %>% 
    mutate(
      p_x = as.double(N_x)/N,
      pmi_xy = p_x*p_y_given_x*log(p_y_given_x/p_y)
    )
  
  tmp5 = tmp4 %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
      # integrate over dy ( well does the trapeziod part of the integration
      # TODO: this integegration could be improved. lead(pmi_xy,default = 0)
      d_I_d_xy = (pmi_xy+lag(pmi_xy,1,default=0))*as.double(y_continuous-lag(y_continuous))/2
    ) %>% group_by(!!!grps,!!!discreteVars) %>% summarise (
      N = max(N),
      N_x = max(N_x),
      p_x = max(p_x),
      I_given_x = sum(d_I_d_xy,na.rm=TRUE)/p_x, # this does the sum over x and computes the integral over y at the same time.
      I_given_x_sd = NA,
      method = "Kernel"
    )
  
  return(tmp5)
    
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - the half window width of the SG filter that smooths the data. This is dependent on data but typically not less that 10.
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_SGolay = function(df, discreteVars, continuousVar, k_05=10, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  
  continuousVar = ensym(continuousVar)
  #joinList = joinList(df, discreteVars)
  
  # N.B. this whole bit of code is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the discrete 
  # NX has group counts (N) and subgroup counts (N_x).
  tmp = df %>% groupwiseCount(discreteVars) %>% labelGroup(discreteVars, x_discrete) %>% mutate(y_continuous=!!continuousVar)
  
  
  
  # if (min(groupSize$N_x) < k_05*2+1) {
  #   return(df %>% group_by(!!!grps) %>% summarise(I = NA))
  # }
  
  # sort by continuousVar
  # assign ungrouped rank (rank)
  # apply grouping discreteVar and sort by continuousVar within group
  # rank within outcome group (groupRank)
  # View(tmp %>% filter(outcome=="low" & test=="K")) # verify collisions get a random rank
  
  
  # get the group wise pdfs (i.e. p_y_given_x) at every observation point in class X
  # process with X grouped and ordered by Y within groups
  # value of Y wrt rank of Y is an inverse CDF function. This is evenly spaced on rank of Y
  # spacing of rank depends on the number of observations (and actually we don't need the value)
  # gradient of this is inverse of PDF. This needs heavy smoothing to be useable subsequently_
  
  
  tmp2 = tmp %>%  group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% group_modify(
    function(d,...) {
      k = k_05*2+1
      samples = min(d$N_x,na.rm=TRUE)
      #if (k >= samples) k = samples-samples%%2-1
      if (k < samples-1) {
        d_xy_d_r = signal::sgolayfilt(d$y_continuous, p=2, n=k, m=1, ts=1.0/samples)
        d_xy_d_r = ifelse(d_xy_d_r <= 0, 0.0001, d_xy_d_r)
        return(
          tibble(
            N_x = d$N_x,
            N = d$N,
            y_continuous = d$y_continuous,
            d_xy_d_r = d_xy_d_r # prevent negative gradient - d_x_d_r is an inverse cdf - always positive
          ) %>% mutate(
            p_y_given_x = 1.0/d_xy_d_r
          )
        )
      } else {
        return(
          tibble(
            N_x = d$N_x,
            N = d$N,
            y_continuous = d$y_continuous,
            d_xy_d_r = rep(NA,length(d$y_continuous)), 
            p_y_given_x = rep(NA,length(d$y_continuous)) 
          )
        )
      }
    }
  )
  
  # get the overall pdf for the combination of all x's (gives us p_y)	at all observation points
  tmp2 = tmp2 %>%  group_by(!!!grps) %>% arrange(y_continuous) %>% group_modify(
    function(d,...) {
      k = k_05*2+1
      samples = min(d$N,na.rm = TRUE)
      # if (k >= samples) k = samples-samples%%2-1
      if (k < samples-1) {
        d_y_d_r = signal::sgolayfilt(d$y_continuous, p=2, n=k, m=1, ts=1.0/samples)
        d_y_d_r = ifelse(d_y_d_r <= 0, 0.0001, d_y_d_r)
        return(
          tibble(
            N_x=d$N_x,
            N = d$N,
            y_continuous = d$y_continuous,
            x_discrete = d$x_discrete,
            p_y_given_x = d$p_y_given_x,
            d_xy_d_r = d$d_xy_d_r,
            d_y_d_r = d_y_d_r # prevent negative gradient - d_x_d_r is an inverse cdf - always positive
          ) %>% mutate(
            p_y = 1.0/d_y_d_r,
            p_x = as.double(N_x)/N
          ) %>% mutate(
            pmi_xy = p_x*p_y_given_x*log(p_y_given_x/p_y)
          )
        )
      } else {
        return(
          tibble(
            N_x=d$N_x,
            N = d$N,
            y_continuous = d$y_continuous,
            x_discrete = d$x_discrete,
            p_y_given_x = d$p_y_given_x,
            d_xy_d_r = d$d_xy_d_r,
            d_y_d_r = rep(NA,length(d$y_continuous)),
            p_y = rep(NA,length(d$y_continuous)),
            p_x = as.double(N_x)/N,
            pmi_xy = rep(NA,length(d$y_continuous)) 
          )
        )
      }
    }
  )
  
  # do the integration
  # points are not evenly spaced in y dimension so piecewise integration
  tmp3 = tmp2 %>% group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% mutate(
    # integrate over dy ( well does the trapeziod part of the integration
    # TODO: this integegration could be improved. lead(pmi_xy,default = 0)
    d_I_d_xy = (pmi_xy+lag(pmi_xy,1,default=0))*as.double(y_continuous-lag(y_continuous))/2
  ) %>% group_by(!!!grps) %>% summarise (
    min_N_x = min(N_x),
    N = max(N),
    min_N_x = min(N_x),
    I = sum(d_I_d_xy,na.rm=TRUE), # this does the sum over x and computes the integral over y at the same time.
    I_sd = NA,
    method = "SGolay"
  ) %>% mutate(
    I = ifelse(min_N_x-1 < k_05*2+1,NA,I)
  ) %>% select(-min_N_x)
  
  
  
  return(tmp3)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is based on the technique described here:
#' 
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' but with the important simplification of using the sliding window K elements wide rather than the k nearest neighbours. This is empirically shown to have little difference on larger datasets
#' and makes this algorithm simple to implement in dbplyr tables.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KWindow = function(df, discreteVars, continuousVar, k_05=4, ...) { 
  grps = df %>% groups()
  
  tmp5 = calculateDiscreteContinuousConditionalMI_KWindow(df, discreteVars, {{continuousVar}}, k_05, ...)
  
  tmp6 = tmp5 %>% group_by(!!!grps) %>% summarise(
    N = max(N),
    I = sum(p_x * I_given_x),
    I_sd = sum(p_x * I_given_x_sd),
    method = "KWindow"
  )
  
  return(tmp6)
}


#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y). I.e. the self information of X conditioned on each of the possible values for X
#' using a sliding window and local entropy measure
#' 
#' This is based on the technique described here:
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' but with the important simplification of using the sliding window K elements wide rather than the k nearest neighbours. This is empirically shown to have little difference on larger datasets
#' and makes this algorithm simple to implement in dbplyr tables.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousConditionalMI_KWindow = function(df, discreteVars, continuousVar, k_05=4L, ...) { 
  k_05 = as.integer(k_05)
  if (k_05<2L) k_05=2L
  grps = df %>% groups()
  continuousVar = ensym(continuousVar)
  joinList = df %>% joinList(discreteVars)
  
  df = df %>% select(!!!grps,!!!discreteVars,!!continuousVar)
  # this is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the categorical X 
  # has group counts (N) and subgroup counts (N_x) 
  
  tmp = df %>% groupwiseCount(discreteVars) %>% mutate(y_continuous=!!continuousVar)
  
  # the knn approach without using neighbours - i.e. a k wide sliding window
  tmp4 = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number())
  tmp4 = tmp4  %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
    
    # correct k for tails of distributions exclusive
    # kRank = row_number(),
    # m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1)+1L,
    # k = lead(kRank,n=k_05,default=max(N_x))-lag(kRank,n=k_05,default=1)+1L
    
    # correct k for tails of distributions inclusive
    # kRank = row_number(),
    # m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1),
    # k = lead(kRank,n=k_05,default=max(N_x))-lag(kRank,n=k_05,default=1)
    
    # dont correct k & exclude tails
    k = local(k_05)*2L,
    m_i = ifelse(N_x < local(k_05)*2L+1L, NA, lead(rank,n=local(k_05))-lag(rank,n=local(k_05)))
    
    # average m_i over 3 window sizes
    # k = k_05*2,
    # m_i = floor((
    #				lead(rank,n=local(k_05))-lag(rank,n=local(k_05))+
    #				lead(rank,n=local(k_05+1L))-lag(rank,n=local(k_05+1L))+
    #				lead(rank,n=local(k_05-1L))-lag(rank,n=local(k_05-1L))
    #				)/3*2)/2
  )  
  
  # there are situations where this estimate can go negative - its to do with the estimation error in small sample sizez
  # and the fact that in the window method the sliding is potentially centred around different points
  # TODO: try and understand this:
  # tmp4 = tmp4 %>% mutate(m_i = ifelse(m_i < k,NA,m_i)) %>% compute()
  tmp4 = tmp4 %>% compute()
  
  # TODO: what happens when there are not enough points in one class? As is we seem to get conditional estimate based on one class
  # but not the other. In the extreme case where there are no points in the minority class we get a degnerate behaviour and
  # an overall mutual information of 0 (which may be technically correct)
  
  
  tmp4 = tmp4 %>% 
    calculateDigamma(k,digammak) %>% 
    calculateDigamma(N,digammaN) %>% 
    calculateDigamma(N_x,digammaN_x) %>% 
    calculateDigamma(m_i,digammam_i) %>% 
    mutate(
      I_i = digammaN-digammaN_x+digammak-digammam_i
    )
  
  tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,!!!discreteVars) %>% summarize(
    N = max(N),
    N_x = max(N_x),
    p_x = max(as.double(N_x)/N),
    I_given_x = mean(I_i,na.rm = TRUE),
    I_given_x_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE))
  ) %>% mutate(method = "KWindow")
  
  
  return(tmp5 %>% group_by(!!!grps))
}










#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a discretisation and discrete MNI estimators
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param bins - the number of bins
#' @param binStrategy - the a function to genereate the bins on a group by group basis
#' @param discreteMethod - What methd will be used to calculate the MI once discretised?
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_DiscretiseByRank = function(df, discreteVars, continuousVar, bins=NA, binStrategy = linearBySize(slope=8,minBins=4,maxBins=100), discreteMethod="Histogram", ...) {
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  
  tmp = df %>% rename(
    y_continuous=!!continuousVar
  )
  
  tmp = tmp %>% discretise_ByRank(y_continuous, y_discrete, bins, binStrategy) %>% compute()
  return(tmp %>% calculateDiscreteDiscreteMI(discreteVars, vars(y_discrete), method=discreteMethod, ...) %>% 
           collect() %>% mutate(method = paste0("DiscretiseByRank - ",method)))
  
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a discretisation and infotheo
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param bins - the number of bins
#' @param binStrategy - the a function to genereate the bins on a group by group basis
#' @param discreteMethod - What methd will be used to calculate the MI once discretised?
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_DiscretiseByValue = function(df, discreteVars, continuousVar, bins=NA, binStrategy = linearBySize(slope=8,minBins=4,maxBins=100), discreteMethod="Histogram", ...) {
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  
  tmp = df %>% rename(
    y_continuous=!!continuousVar
  )
  
  tmp = tmp %>% discretise_ByValue(y_continuous, y_discrete, bins, binStrategy) %>% compute()
  
  return(tmp %>% calculateDiscreteDiscreteMI(discreteVars, vars(y_discrete), method=discreteMethod, ...) %>% 
           collect() %>% mutate(method = paste0("DiscretiseByValue - ",method)))
  
  #return(
  #  df %>% group_modify(function(d,...) {
  #      d = d %>% mutate(x_discrete = as.integer(as.factor(!!discreteVar))) %>% select(x_discrete, y_continuous = !!continuousVar) %>% infotheo::discretize(n=bins)
  #     return(data.frame(I=infotheo::mutinformation(d$x_discrete, d$y_continuous)))
  # })
  #)
}


#### ----

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is an implementation of the technique described here:
#' 
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' 
#' But it is very slow. Empirically it also does not give any better estimate that the KWindow method.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @param useKWindow - will switch to using the much faster KWindow estimator for larger sample sizes (>500) when the difference between the 2 methods is negligable
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KNN = function(df, discreteVars, continuousVar, k_05=4L, useKWindow = TRUE,...) {
  grps = df %>% groups()

  maxN = df %>% group_by(!!!grps) %>% summarise(N = n()) %>% summarise(maxN = max(N, na.rm=TRUE)) %>% pull(maxN)
  
  if(maxN > 500 && useKWindow) {
    message("using KWindow method instead of KNN for larger sample")
    return(calculateDiscreteContinuousMI_KWindow(df, discreteVars, {{continuousVar}}, k_05))
  }
  
  tmp5 = calculateDiscreteContinuousConditionalMI_KNN(df, discreteVars, {{continuousVar}}, k_05, useKWindow, ...)

  tmp6 = tmp5 %>% group_by(!!!grps) %>% summarise(
    N = max(N),
    I = sum(p_x * I_given_x),
    I_sd = sum(p_x * I_given_x_sd),
    method = "KNN"
  )
  
  return(tmp6)
}


#' calculate pointwise mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is an implementation of the technique described here:
#' 
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' 
#' But it is very slow. Empirically it also does not give any better estimate that the KWindow method.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @param useKWindow - will switch to using the much faster KWindow estimator for larger sample sizes (>500) when the difference between the 2 methods is negligable
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousConditionalMI_KNN = function(df, discreteVars, continuousVar, k_05=4L, useKWindow = TRUE,...) {
  grps = df %>% groups()
  maxN = df %>% group_by(!!!grps) %>% summarise(N = n()) %>% summarise(maxN = max(N, na.rm=TRUE)) %>% pull(maxN)
  if(maxN > 500 && useKWindow) {
    message("using KWindow method instead of KNN for larger sample")
    return(calculateDiscreteContinuousConditionalMI_KWindow(df, discreteVars, {{continuousVar}}, k_05))
  }
  k_05 = as.integer(k_05)
  if (k_05<2L) k_05=2L
  continuousVar = ensym(continuousVar)
  
  joinList = df %>% joinList(discreteVars)
  
  df = df %>% select(!!!grps,!!!discreteVars,!!continuousVar)
  # this is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the categorical X 
  # NX has group counts (N) and subgroup counts (N_x) 
  tmp = df %>% groupwiseCount(discreteVars) %>% labelGroup(discreteVars, x_discrete) %>% mutate(y_continuous=!!continuousVar)
  
  # add row numbers of sequence - sorted by value, and sorted by value in groups
  tmp = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number()) %>% 
    group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(groupRank = row_number()) %>% compute()
  
  tmp = tmp %>% group_by(!!!grps,!!!discreteVars) %>% arrange(y_continuous) %>% mutate(
    rankMax = lead(rank,n=k_05*2L),
    rankMin = lag(rank,n=k_05*2L,1L)
  ) %>% mutate(rankMax = ifelse(is.na(rankMax),N,rankMax))
  
  # list of join variables for join by value
  join2List = df %>% joinList(defaultJoin="join")
  
  # rhs of join by value
  tmp_join = tmp %>% ungroup() %>% rename(
    y_continuous_knn = y_continuous,
    x_discrete_knn = x_discrete,
    rank_knn = rank) %>%
    select(!!!grps,x_discrete_knn,y_continuous_knn,rank_knn) %>% mutate(join=1L) %>% compute()
  
  # TODO: this is unuseably inefficient
  
  tmp4 = tmp %>% mutate(join=1) %>% inner_join(tmp_join, by=join2List) %>% 
    filter(rankMax >= rank_knn & rankMin <= rank_knn) %>%
    mutate(y_diff = abs(y_continuous - y_continuous_knn)) %>% 
    group_by(!!!grps,rank) %>% 
    arrange(y_diff) %>% 
    mutate(sameGroup=ifelse(x_discrete_knn==x_discrete,1L,0L)) %>% #, differentGroup=ifelse(x_discrete_knn==x_discrete,0,1)) %>%
    mutate(kDist = cumsum(sameGroup), m_i = row_number()) %>%
    filter(kDist == local(k_05*2L+1L) & sameGroup==1L) %>% 
    mutate(
      k = ifelse(N < k_05*2L+1L, NA, kDist-1L), 
      m_i = m_i-1
    ) %>% compute() # this is the strange definintion of knn in the original paper
  
  tmp4 = tmp4 %>% 
    calculateDigamma(k,digammak) %>% 
    calculateDigamma(N,digammaN) %>% 
    calculateDigamma(N_x,digammaN_x) %>% 
    calculateDigamma(m_i,digammam_i) %>% 
    mutate(
      I_i = digammaN-digammaN_x+digammak-digammam_i
    )
  
  tmp5 = tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps,!!!discreteVars) %>% summarize(
    N = max(N),
    N_x = max(N_x),
    p_x = max(as.double(N_x)/N),
    I_given_x = mean(I_i,na.rm = TRUE),
    I_given_x_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE))
  ) %>% mutate(method = "KNN")
  
  return(tmp5)
}




#### ----

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a compression algorithm based on.
#' 
#' Universal and accessible entropy estimation using a compression algorithm Ram Avinery, Micha Kornreich, Roy Beck https://arxiv.org/abs/1709.10164
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_Compression = function(df, discreteVars, continuousVar, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  continuousVar = ensym(continuousVar)
  
  tmp = df %>% rename(
    y_continuous=!!continuousVar
  ) 
  
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  # list of join variables for join by value
  join2List = c(grpsList, "join")
  
  # define the discretisation strategy from 4 to 256 bins depending on sample size
  # Cy defines the maximum entropy possible in that discretised 
  # map continuous variable to a stream of bytes
  
  tmp2 = tmp %>% group_by(!!!grps) %>% discretise_ByRank(y_continuous,y_discrete,binStrategy = linearBySize(8,4,256)) %>% mutate(
    y_raw = as.raw(y_discrete-1)
  )
  
  # https://www.sqlshack.com/compression-decompression-functions-sql-server-2016/
  # private function
  compressionEntropy = function(d,g,...) {
    Cy = n_distinct(d$y_discrete)
    #classes = unique(d$y_discrete)
    N = length(d$y_raw)
    C0 = length(memCompress(as.raw(rep(0,N))))
    # https://math.stackexchange.com/questions/1406747/whats-theoretical-maximum-information-compression-rate
    # in this case our we are looking to find L (compressed message length) based on max possible entropy for message type (log(Cy)) and uncompressed message length (N)
    C1 = as.double(C0+N)/log(Cy) #max(sapply(c(1:10),function(i) length(memCompress(as.raw(sample(classes,size=N,replace=TRUE)-1)))))
    C = length(memCompress(as.vector(d$y_raw)))
    if (C > C1) C=C1 # prevent entropy exceeding theoretical maximum
    H = as.double(C-C0)/(C1-C0) * log(Cy) # original paper includes a degrees of freedom parameter here. with this setup this can only be one...?
    #TODO: the log Cy term here is implicated in differences between this and other algorithms
    return(tibble(
      H=H,
      N=N
    ))
  }
  
  # overall entropy of Y
  tmp3 = tmp2 %>% group_by(!!!grps) %>% group_modify(compressionEntropy) %>% rename(Hy = H) %>% mutate(join = 1)
    # Hy = calculateEntropy(y_raw, method=method),
    # Hy = compressionEntropy(y_raw),
    # Cy = n_distinct(y_raw)),
    # Cy = max(as.integer(y_raw))-min(as.integer(y_raw))+1,
    # Cy = 256, # y is continuous therefore can have any number of classes
    # N = n()
  # ) %>% mutate(join=1)
  
  # entropy of Y given X=x
  tmp4 = tmp2 %>% group_by(!!!grps,!!!discreteVars) %>% group_modify(compressionEntropy) %>% rename(NX = N, Hygivenx = H) %>% 
    group_by(!!!grps) %>% mutate(N = sum(NX))
  
    # NX = n(),
    # Hygivenx = compressionEntropy(y_raw)
    # Hygivenx = calculateEntropy(y_raw,Cy,method)
  # ) %>% mutate(join=1)
  
  tmp5 =  tmp4 %>% summarise(
      Hygivenx = sum( as.double(NX) / N * Hygivenx )
    ) %>% mutate(join = 1)
  
  # https://en.wikipedia.org/wiki/Conditional_entropy
  # H(Y|X) = SUM(over all x) P(x)*H(Y|X=x)
  # I(X,Y) = H(Y) - H(Y|X)
  tmp6 = tmp3 %>% inner_join(tmp5, by=join2List) %>% ungroup() %>% group_by(!!!grps,N)  %>% summarise(
    I = Hy - Hygivenx,
    I_sd = as.numeric(NA),
    method = "Compression"
  )
  
  return(tmp6)
  
}

#' calculate conditional mutual information between a discrete value (X) and a continuous value (Y) using estimates of differential entropy
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @param entropyMethod - the method used to calculate the entropy (see ?tidyinfostats::calculateEntropy) - defaults to "Grassberger"
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteContinuousConditionalMI_Entropy = function(df, discreteVars, continuousVar, entropyMethod="Quantile", ...) {
  stop("not working")
  grps = df %>% groups()
  #joinList = df %>% joinList(defaultJoin = "join")
  # list of join variables for join by value
  innerJoinList = df %>% joinList(discreteVars)
  outerJoinList = df %>% joinList(defaultJoin = "tmp_join")
  continuousVar = ensym(continuousVar)
  
  tmp_H_y = df %>% group_by(!!!grps) %>% 
    calculateContinuousEntropy(!!continuousVar, method = entropyMethod, ...) %>% 
    rename(I_y = I, I_y_sd = I_sd) %>% mutate(tmp_join = 1L)
  
  tmp = df %>% 
    group_by(!!!grps, !!!discreteVars) %>% 
    calculateContinuousEntropy(!!continuousVar, method = entropyMethod, ...) %>% 
    rename(I_given_x = I, I_given_x_sd = I_sd)
  
  tmp_p_x = df %>% group_by(!!!grps) %>% 
    groupwiseCount(discreteVars, summarise = TRUE) %>% 
    mutate(p_x=as.double(N_x)/N)
  
  tmp2 = tmp %>% left_join(tmp_p_x, by=innerJoinList) %>% mutate(tmp_join = 1L)
  tmp2 = tmp2 %>% left_join(tmp_H_y, by=outerJoinList)
  
  tmp5 = tmp2 %>% 
    group_by(!!!grps,!!!discreteVars) %>% 
    mutate(
      I_given_x = I_y*p_x - I_given_x,
      I_given_x_sd = I_y_sd*p_x + I_given_x_sd,
      method = paste0("Entropy - ",entropyMethod)
    )
  browser()
  return(tmp5)
  
}


#' calculate mutual information between a discrete value (X) and a continuous value (Y) using estimates of differential entropy
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param discreteVars - the column(s) of the categorical value (X) quoted by vars(...)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @param entropyMethod - the method used to calculate the entropy (see ?tidyinfostats::calculateEntropy) - defaults to "Grassberger"
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteContinuousMI_Entropy = function(df, discreteVars, continuousVar, entropyMethod="Quantile", ...) {
  grps = df %>% groups()
  return(
    df %>% calculateDiscreteContinuousConditionalMI_Entropy(discreteVars, {{continuousVar}}, entropyMethod, ...) %>% 
      group_by(!!!grps) %>%
      summarise(
        I = sum(p_x*I_given_x),
        I_sd = sum(I_given_x_sd*p_x)
      )
  )
}


