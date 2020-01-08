#' calculate mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param groupXVar - the column of the categorical value (X)
#' @param valueYVar - the column of the continuous value (Y)
#' @param method - the method employed - valid options are "KWindow","KNN","SGolay","DiscreteByRank","DiscreteByValue","Compression","Entropy".
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI = function(df, groupXVar, valueYVar, method="KWindow", ...) {
  switch (method,
    KWindow = calculateDiscreteContinuousMI_KWindow(df, {{groupXVar}}, {{valueYVar}}, ...),
    KNN = calculateDiscreteContinuousMI_KNN(df, {{groupXVar}}, {{valueYVar}}, ...),
    SGolay = calculateDiscreteContinuousMI_SGolay(df, {{groupXVar}}, {{valueYVar}}, ...),
    DiscretiseByRank = calculateDiscreteContinuousMI_DiscretiseByRank(df, {{groupXVar}}, {{valueYVar}}, ...),
    DiscretiseByValue = calculateDiscreteContinuousMI_DiscretiseByValue(df, {{groupXVar}}, {{valueYVar}}, ...),
    Compression = calculateDiscreteContinuousMI_Compression(df, {{groupXVar}}, {{valueYVar}}, ...),
    {stop(paste0(method," not a valid option"))}
    # Entropy = calculateDiscreteContinuousMI_Entropy(df, {{groupXVar}}, {{valueYVar}}, method="Empirical", ...),
  )
}

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

#' calculate mutual information between a categorical value (X) and a continuous value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - the half window width of the SG filter that smooths the data. This is dependent on data but typically not less that 10.
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_SGolay = function(df, discreteVar, continuousVar, k_05=10, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)

  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList, as.character(discreteVar))
  
  # N.B. this whole bit of code is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the discrete 
  # NX has group counts (N) and subgroup counts (N_X).
  groupSize = df %>% group_by(!!!grps,!!discreteVar) %>% summarise(N_X = n()) %>% group_by(!!!grps) %>% mutate(N = sum(N_X))
  
  # if (min(groupSize$N_X) < k_05*2+1) {
  #   return(df %>% group_by(!!!grps) %>% summarise(I = NA))
  # }
  
  # sort by valueYVar
  # assign ungrouped rank (rank)
  # apply grouping groupXVar and sort by valueYVar within group
  # rank within outcome group (groupRank)
  # View(tmp %>% filter(outcome=="low" & test=="K")) # verify collisions get a random rank
  
  tmp = df %>% inner_join(groupSize, by=joinList) %>% mutate(
    y_continuous=!!continuousVar,
    x_discrete=!!discreteVar)
  # get the group wise pdfs (i.e. p_y_given_x) at every observation point in class X
  # process with X grouped and ordered by Y within groups
  # value of Y wrt rank of Y is an inverse CDF function. This is evenly spaced on rank of Y
  # spacing of rank depends on the number of observations (and actually we don't need the value)
  # gradient of this is inverse of PDF. This needs heavy smoothing to be useable subsequently_
  
  
  tmp2 = tmp %>%  group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% group_modify(
    function(d,...) {
      k = k_05*2+1
      samples = min(d$N_X)
      #if (k >= samples) k = samples-samples%%2-1
      if (k < samples-1) {
        return(
          tibble(
            N_X = d$N_X,
            N = d$N,
            y_continuous = d$y_continuous,
            d_xy_d_r = signal::sgolayfilt(d$y_continuous, p=2, n=k, m=1, ts=1.0/samples)
          ) %>% mutate(
            p_y_given_x = 1.0/d_xy_d_r
          )
        )
      } else {
        return(
          tibble(
            N_X = d$N_X,
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
      samples = min(d$N)
      # if (k >= samples) k = samples-samples%%2-1
      if (k < samples-1) {
        return(
          tibble(
            N_X=d$N_X,
            N = d$N,
            x_discrete = d$x_discrete,
            y_continuous = d$y_continuous,
            p_y_given_x = d$p_y_given_x,
            d_xy_d_r = d$d_xy_d_r,
            d_y_d_r = signal::sgolayfilt(d$y_continuous, p=2, n=k, m=1, ts=1.0/samples)
          ) %>% mutate(
            p_y = 1.0/d_y_d_r,
            p_x = as.double(N_X)/N
          ) %>% mutate(
            pmi_xy = p_x*p_y_given_x*log(p_y_given_x/p_y)
          )
        )
      } else {
        return(
          tibble(
            N_X=d$N_X,
            N = d$N,
            x_discrete = d$x_discrete,
            y_continuous = d$y_continuous,
            p_y_given_x = d$p_y_given_x,
            d_xy_d_r = d$d_xy_d_r,
            d_y_d_r = rep(NA,length(d$y_continuous)),
            p_y = rep(NA,length(d$y_continuous)),
            p_x = as.double(N_X)/N,
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
    d_I_d_xy = (pmi_xy+lag(pmi_xy,1,default=0))*as.double(y_continuous-lag(y_continuous))/2
  ) %>% group_by(!!!grps) %>% summarise (
    min_N_X = min(N_X),
    I = sum(d_I_d_xy,na.rm=TRUE),
    I_sd = NA,
    method = "SGolay"
  ) %>% mutate(
    I = ifelse(min_N_X-1 < k_05*2+1,NA,I)
  ) %>% select(-min_N_X)
  
  
  
  return(tmp3)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a sliding window and local entropy measure
#' 
#' This is based on the technique described here:
#' B. C. Ross, “Mutual information between discrete and continuous data sets,” PLoS One, vol. 9, no. 2, p. e87357, Feb. 2014 [Online]. Available: http://dx_doi.org/10.1371/journal.pone.0087357
#' but with the important simplification of using the sliding window K elements wide rather than the k nearest neighbours. This is empirically shown to have little difference on larger datasets
#' and makes this algorithm simple to implement in dbplyr tables.
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KWindow = function(df, discreteVar, continuousVar, k_05=4, ...) { 
  k_05 = as.integer(k_05)
  if (k_05<2) k_05=2
  grps = df %>% groups()
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)
  
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList, as.character(discreteVar))
  
  df = df %>% select(!!!grps,!!discreteVar,!!continuousVar)
  # this is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the categorical X 
  # NX has group counts (N) and subgroup counts (N_X) 
  grpCounts = df %>% group_by(!!!grps,!!discreteVar) %>% summarise(N_X = n()) %>% group_by(!!!grps) %>% mutate(N = sum(N_X, na.rm=TRUE)) %>% compute()
  
  tmp = df %>% inner_join(grpCounts, by=joinList) %>% mutate(
    x_discrete=!!discreteVar,
    y_continuous=!!continuousVar)
  
  # the knn approach without using neighbours - i.e. a k wide sliding window
  tmp4 = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number())
  tmp4 = tmp4 %>% group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% mutate(
    
    # correct k for tails of distributions exclusive
    # kRank = row_number(),
    # m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1)+1L,
    # k = lead(kRank,n=k_05,default=max(N_X))-lag(kRank,n=k_05,default=1)+1L
    
    # correct k for tails of distributions inclusive
    # kRank = row_number(),
    # m_i = lead(rank,n=k_05,default=max(N))-lag(rank,n=k_05,default=1),
    # k = lead(kRank,n=k_05,default=max(N_X))-lag(kRank,n=k_05,default=1)
    
    # dont correct k & exclude tails
    k = k_05*2,
    m_i = lead(rank,n=k_05)-lag(rank,n=k_05)
    
    # average m_i over 3 window sizes
    # k = k_05*2,
    # m_i = floor((
	  #				lead(rank,n=local(k_05))-lag(rank,n=local(k_05))+
	  #				lead(rank,n=local(k_05+1L))-lag(rank,n=local(k_05+1L))+
	  #				lead(rank,n=local(k_05-1L))-lag(rank,n=local(k_05-1L))
	  #				)/3*2)/2
  )  %>% compute()
  
  tmp4 = tmp4 %>% 
    calculateDigamma(k,digammak) %>% 
    calculateDigamma(N,digammaN) %>% 
    calculateDigamma(N_X,digammaN_X) %>% 
    calculateDigamma(m_i,digammam_i) %>% 
    mutate(
      I_i = digammaN-digammaN_X+digammak-digammam_i
    )
  
  tmp5 = tryCatch(
    tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps) %>% summarize(
      # TODO: catch the NaN when no group is larger that k_05*2+1
      I = mean(I_i,na.rm = TRUE),
      I_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE)),
      method = "KWindow"
    ),
    warning = function(w) {
      tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps) %>% summarize(
        I = NA,
        I_sd = NA,
        method = "KWindow"
      )
    }
  )
  return(tmp5)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a discretisation and discrete MNI estimators
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param bins - the number of bins
#' @param binStrategy - the a function to genereate the bins on a group by group basis
#' @param discreteMethod - What methd will be used to calculate the MI once discretised?
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_DiscretiseByRank = function(df, discreteVar, continuousVar, bins=NA, binStrategy = linearBySize(slope=8,minBins=4,maxBins=100), discreteMethod="Histogram", ...) {
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  
  tmp = df %>% rename(
    x_discrete=!!discreteVar,
    y_continuous=!!continuousVar
  )
  
  tmp = tmp %>% discretise_ByRank(y_continuous, y_discrete, bins, binStrategy) %>% compute()
  return(tmp %>% calculateDiscreteDiscreteMI(x_discrete, y_discrete, method=discreteMethod, ...) %>% collect() %>% mutate(method = paste0("DiscretiseByRank - ",method)))
  
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a discretisation and infotheo
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param bins - the number of bins
#' @param binStrategy - the a function to genereate the bins on a group by group basis
#' @param discreteMethod - What methd will be used to calculate the MI once discretised?
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_DiscretiseByValue = function(df, discreteVar, continuousVar, bins=NA, binStrategy = linearBySize(slope=8,minBins=4,maxBins=100), discreteMethod="Histogram", ...) {
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)
  grps = df %>% groups()
  
  tmp = df %>% rename(
    x_discrete=!!discreteVar,
    y_continuous=!!continuousVar
  )
  
  tmp = tmp %>% discretise_ByValue(y_continuous, y_discrete, bins, binStrategy) %>% compute()
  
  return(tmp %>% calculateDiscreteDiscreteMI(x_discrete, y_discrete, method=discreteMethod, ...) %>% collect() %>% mutate(method = paste0("DiscretiseByValue - ",method)))
  
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
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param k_05 - half the sliding window width - this should be a small number like 1,2,3.
#' @param useKWindow - will switch to using the much faster KWindow estimator for larger sample sizes (>500) when the difference between the 2 methods is negligable
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_KNN = function(df, discreteVar, continuousVar, k_05=4, useKWindow = TRUE,...) { #a=0.992, b=1) {
  grps = df %>% groups()
  maxN = df %>% group_by(!!!grps) %>% summarise(N = n()) %>% summarise(maxN = max(N, na.rm=TRUE)) %>% pull(maxN)
  if(maxN > 500 && useKWindow) {
    # warning("using KWindow method instead of KNN for larger sample")
    return(calculateDiscreteContinuousMI_KWindow(df,{{discreteVar}}, {{continuousVar}}, k_05))
  }
  #if (length(df) > 500) df = sample(df,size=500)
  k_05 = as.integer(k_05)
  if (k_05<2) k_05=2
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)
  
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  joinList = c(grpsList, as.character(discreteVar))
  
  df = df %>% select(!!!grps,!!discreteVar,!!continuousVar)
  # this is confusing because groups mean 2 things here - the 
  # different types of Y (grps) which should be preserved and the categorical X 
  # NX has group counts (N) and subgroup counts (N_X) 
  grpCounts = df %>% group_by(!!!grps,!!discreteVar) %>% summarise(N_X = n()) %>% group_by(!!!grps) %>% mutate(N = sum(N_X)) %>% compute()
  
  # add in overall counts of groups
  tmp = df %>% inner_join(grpCounts, by=joinList) %>% mutate(
    x_discrete=!!discreteVar,
    y_continuous=!!continuousVar)
  
  # add row numbers of sequence - sorted by value, and sorted by value in groups
  tmp = tmp %>% group_by(!!!grps) %>% arrange(y_continuous) %>% mutate(rank = row_number()) %>% 
    group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% mutate(groupRank = row_number()) %>% compute()
  
  tmp = tmp %>% group_by(!!!grps,x_discrete) %>% arrange(y_continuous) %>% mutate(
    rankMax = lead(rank,n=k_05*2),
    rankMin = lag(rank,n=k_05*2,1)
  ) %>% mutate(rankMax = ifelse(is.na(rankMax),N,rankMax))
  
  # list of join variables for join by value
  join2List = c(grpsList, "join")
  
  # rhs of join by value
  tmp_join = tmp %>% mutate(join=1) %>% rename(
    y_continuous_knn = y_continuous, 
    x_discrete_knn = x_discrete,
    rank_knn = rank) %>%
    select(!!!grps,join,y_continuous_knn,x_discrete_knn, rank_knn) %>% compute()
  
  # TODO: this is unuseably inefficient
  
  tmp4 = tmp %>% mutate(join=1) %>% inner_join(tmp_join, by=join2List) %>% 
    filter(rankMax >= rank_knn & rankMin <= rank_knn) %>%
    mutate(y_diff = abs(y_continuous - y_continuous_knn)) %>% 
    group_by(!!!grps,rank) %>% 
    arrange(y_diff) %>% 
    mutate(sameGroup=ifelse(x_discrete_knn==x_discrete,1L,0L)) %>% #, differentGroup=ifelse(x_discrete_knn==x_discrete,0,1)) %>%
    mutate(kDist = cumsum(sameGroup), m_i = row_number()) %>%
    filter(kDist == local(k_05*2L+1) & sameGroup==1) %>% 
    mutate(k = kDist-1, m_i=m_i-1) %>% compute() # this is the strange definintion of knn in the original paper
  
  tmp4 = tmp4 %>% 
    calculateDigamma(k,digammak) %>% 
    calculateDigamma(N,digammaN) %>% 
    calculateDigamma(N_X,digammaN_X) %>% 
    calculateDigamma(m_i,digammam_i) %>% 
    mutate(
      I_i = digammaN-digammaN_X+digammak-digammam_i
    )
  
  tmp5 = tryCatch(tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps) %>% summarize(
    I = mean(I_i,na.rm = TRUE),
    I_sd = sd(I_i,na.rm = TRUE)/sqrt(max(N,na.rm=TRUE)),
    method = "KNN"
  ),
  warning = function(w) {
    tmp4 %>% filter(!is.na(I_i)) %>% group_by(!!!grps) %>% summarize(
      I = NA,
      I_sd = NA,
      method = "KNN"
    )
  }
  )
  return(tmp5)
}




#### ----

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a compression algorithm based on.
#' 
#' Universal and accessible entropy estimation using a compression algorithm Ram Avinery, Micha Kornreich, Roy Beck https://arxiv.org/abs/1709.10164
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param discreteVar - the column of the categorical value (X)
#' @param continuousVar - the column of the continuous value (Y)
#' @param collect - unless TRUE this function will fail on dbplyr tables as there is no SQL implementation
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteContinuousMI_Compression = function(df, discreteVar, continuousVar, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  discreteVar = ensym(discreteVar)
  continuousVar = ensym(continuousVar)
  
  tmp = df %>% rename(
    y_continuous=!!continuousVar,
    x_discrete=!!discreteVar
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
  tmp4 = tmp2 %>% group_by(!!!grps,x_discrete) %>% group_modify(compressionEntropy) %>% rename(NX = N, Hygivenx = H) %>% 
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
  tmp6 = tmp3 %>% inner_join(tmp5, by=join2List) %>% ungroup() %>% group_by(!!!grps)  %>% summarise(
    I = Hy - Hygivenx,
    I_sd = NA,
    method = "Compression"
  )
  
  return(tmp6)
  
}


