#' calculate entropy of a sequence of 
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param method - the method employed - valid options are "MontgomerySmith", "Histogram", "Grassberger", "InfoTheo", "Compression"
#' @param ... - the other parameters are passed onto the implementations
#' @return a single value for the entropy of the vector
#' @import dplyr
#' @export
calculateDiscreteEntropy = function(df, groupVars, method, ...) {
  switch (method,
          MontgomerySmith = calculateDiscreteEntropy_MontgomerySmith(df, groupVars, ...),
          Histogram = calculateDiscreteEntropy_Histogram(df, groupVars, ...),
          Grassberger = calculateDiscreteEntropy_Grassberger(df, groupVars, ...),
          InfoTheo = calculateDiscreteEntropy_InfoTheo(df, groupVars, ...),
          Compression = calculateDiscreteEntropy_Compression(df, groupVars, ...),
          {stop(paste0(method," not a valid option"))}
  )
}



#' calculate entropy of an optionally ordered discrete value (X) using estimates of entropy from method 2 in the following:
#' 
#' S. Montgomery-Smith and T. Schürmann, “Unbiased Estimators for Entropy and Class Number,” arXiv, 18-Oct-2014. Available: http://arxiv.org/abs/1410.5002
#' 
#' This is not particularly useful as requires large sample size before it becomes accurate.
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param orderingVar - (optional) the column of an ordering variable (e.g. time) - if missing assumes df order,
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteEntropy_MontgomerySmith = function(df, groupVars, orderingVar = NULL, ...) {
  
  # Euler–Mascheroni constant (lambda)
  lambda = 0.577215664901532
  
  grps = df %>% groups()
  # groupVars = ensyms(groupVars)
  
  if (is.null(orderingVar)) {
    # seq here is to fool dbplyr queries into a valid sql syntax as row_number window queries are non-deterministic unless there is some intrinsic ordering.
    tmp = df %>% mutate(seq = 1)
  } else {
    orderingVar = ensym(orderingVar)
    tmp = df %>% mutate(seq = !!orderingVar)
  }
  
  # define rank and top level counts
  tmp = tmp %>% group_by(!!!grps) %>% arrange(seq) %>% mutate(
    rank = row_number(),
    N = n()
  ) %>% groupMutate(
    C_x = n_distinct(!!!groupVars)
  )
  
  # define groupwise rank delta and count
  # tmp = tmp %>% group_by(!!!grps, !!!groupVars) %>% arrange(seq) %>% mutate( # would also probably work
  # browser()
  tmp = tmp %>% group_by(!!!grps, !!!groupVars) %>% arrange(rank) %>% mutate(
    k = lead(rank,1,NA)-rank, # this is the N in the paper - j is as mentioned in the paper,
    # TODO: they discuss a esimator be redone for different values of lead and averaged.
    # I don't understand how this doesn't just increase the estimate as the large j is the larger digamma k is.
    #k2 = lead(rank,2,NA)-rank, would need digamma calc & average
    #k3 = lead(rank,3,NA)-rank,
    NX = n()
  ) %>% mutate(k = ifelse(is.na(k), NX-rank+1, k))
  # if k is not known assume that it is the next one. This "fills in" unknown values
  
  # digamma(n) = Harmonic(n-1)-lambda
  # Entropy estimated by Harmonic(n-1) = digamma(n)+lambda
  
  tmp = tmp %>% calculateDigamma(k, digammak) %>% 
    #calculateDigamma(NX, digammaNX) %>% 
    #calculateDigamma(N, digammaN) %>% 
    #calculateDigamma(C_x, digammaC_x) %>% 
    #mutate(Hx = digammaN-digammaNX+lambda-digammak) #TODO: some adjustment here is possivle
    #mutate(Hx = -digammaN+digammaNX+lambda+digammak) #TODO: some adjustment here is possivle - this is zero based
    mutate(Hx = digammak+lambda)
  
  tmp2 = tmp %>% group_by(!!!grps, N) %>% summarise(
    # H = (mean(digammak, na.rm=TRUE) + lambda), #-digammaC_x+log(C_x), na.rm=TRUE) + lambda),
    I = mean(Hx, na.rm=TRUE),
    # H_sd = sd(digammak, na.rm=TRUE)/sqrt(max(N, na.rm=TRUE)) #-digammaC_x+log(C_x), na.rm=TRUE)
    I_sd = sd(Hx, na.rm=TRUE)/sqrt(N), #-digammaC_x+log(C_x), na.rm=TRUE)
    method = "MontgomerySmith" 
  )  
  
  # %>% mutate(
  #  H = H/log(C_x),
  #  H_sd = H_sd/log(C_x) #TODO: making this up
  #)
  # browser()
  return(tmp2 %>% compute())
  
}



#' calculate entropy of an optionally discrete value (X) using a histogram approach
#' 
#' see: R. de Matos Simoes and F. Emmert-Streib, “Influence of statistical estimators of mutual information and data heterogeneity on the inference of gene regulatory networks,” PLoS One, vol. 6, no. 12, p. e29279, Dec. 2011 [Online]. Available: http://dx.doi.org/10.1371/journal.pone.0029279
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @param mm - Apply a miller-madow adjustment to the result
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteEntropy_Histogram = function(df, groupVars,  countVar=NULL, mm=TRUE, ...) {
  grps = df %>% groups()
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  # groupVars = ensyms(groupVars)
  
  tmp = df %>% calculateSelfInformation_Histogram(groupVars, !!countVar, mm)
  
  tmp3 = tmp %>% ungroup() %>% group_by(!!!grps, N) %>% summarise(
    I = sum(p_x*I_x,na.rm = TRUE), 
    I_sd = NA,
    method = "Histogram"
  )
  
  return(tmp3 %>% ungroup())
}

#' calculate entropy of an optionally discrete value (X) using a histogram approach using the following method
#' 
#' P. Grassberger, “Entropy Estimates from Insufficient Samplings,” arXiv [physics.data-an], 29-Jul-2003 [Online]. Available: http://arxiv.org/abs/physics/0307138
#' 
#' but with a digamma based function (rather than harmonics) detailed in eqns 31 & 35.
#' For our purposes we fix l=0 to give the form in eqn 27. The error in this method is supposedly better for undersampled cases (where number of bins similar to number of samples)
#'  
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteEntropy_Grassberger = function(df, groupVars, countVar=NULL, ...) {
  grps = df %>% groups()
  countVar = tryCatch(ensym(countVar),error = function(e) NULL)
  # groupVars = ensyms(groupVars)
  
  tmp = df %>% groupwiseCount(groupVars, !!countVar, summarise=TRUE) %>% group_by(!!!grps) %>% mutate(C_x = n())
  tmp2 = tmp %>% calculateDigamma(N_x,digamma_N_x) %>% 
    mutate(
      p_x = as.double(N_x)/N,
      G_N_x = digamma_N_x + ((-1)^N_x)/(N_x*(N_x+1)), #TODO: consider expanding for l=1, l=2, etc...
    )
  
  tmp3 = tmp2 %>% ungroup() %>% group_by(!!!grps, N) %>% summarise(
    I = log(N) - sum(p_x * G_N_x,na.rm = TRUE),
    I_sd = NA,
    method = "Grassberger"
  ) 
  
  return(tmp3 %>% ungroup())
}




#' calculate entropy of an optionally discrete value (X) using a infotheo library
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the column of the discrete value (X)
#' @param infoTheoMethod - the method of the entropy estimator ("mm","emp","shrink","sg")
#' @param collect - if false (the default) this will fail on a dbplyr table as this is not supported in SQL
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteEntropy_InfoTheo = function(df, groupVars, infoTheoMethod="mm", collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  # groupVars = ensyms(groupVars)
  
  groupsJoin = as.vector(sapply(groupVars,as_label))
  groupIds = df %>% ungroup() %>% select(!!!groupVars) %>% distinct() %>% arrange(!!!groupVars) %>% mutate(
    # raw_x = TODO convert groupVars to integer... Doesn't matter how - ignores groups even, although maybe shouldn't
    x_int = row_number()
  )
  
  tmp = df %>% ungroup() %>% group_by(!!!grps) %>% #mutate(
    #N = n(),
    #C_x = n_distinct(!!!groupVars),
    # raw_x = TODO convert groupVars to integer... Doesn't matter how - ignores groups even, although maybe shouldn't
  #) %>% 
    left_join(groupIds, by = groupsJoin)
  
  tmp2 = tmp %>% ungroup() %>% group_by(!!!grps,N) %>% group_modify(function(d,...) {
    tibble(
      I = infotheo::entropy(d$x_int, method=infoTheoMethod),
      I_sd = NA,
      method = "InfoTheo"
    )
  })
  
  return(tmp2)
}

#' calculate mutual information between a categorical value (X) and a continuous value (Y) using a compression algorithm based on.
#' 
#' Universal and accessible entropy estimation using a compression algorithm Ram Avinery, Micha Kornreich, Roy Beck https://arxiv.org/abs/1709.10164
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupVars - the column of the discrete value (X)
#' @param orderingVar - (optional) the column of an ordering variable (e.g. time) - if missing assumes df order,
#' @param collect - if TRUE will collect dbplyr tables before processing, otherwise (the default) will fail on dbplyr tables
#' @return a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
#' @import dplyr
#' @export
calculateDiscreteEntropy_Compression = function(df, groupVars, orderingVar = NULL, collect=FALSE, ...) {
  df = collectDf(df,collect)
  grps = df %>% groups()
  # groupVars = ensyms(groupVars)
  if (length(grps)==0) {
    grpsList = NULL
  } else {
    grpsList = sapply(grps,as.character)
  }
  groupsJoin = c(grpsList,as.vector(sapply(groupVars,as_label)))
  
  if (is.null(orderingVar)) {
    # seq here is to fool dbplyr queries into a valid sql syntax as row_number window queries are non-deterministic unless there is some intrinsic ordering.
    tmp = df %>% mutate(tmp_seq=1)
  } else {
    orderingVar = ensym(orderingVar)
    tmp = df %>% mutate(tmp_seq=!!orderingVar)
  }
  
  # convert discrete values of x (defined as combination of groupVars) to integers.
  groupIds = df %>% ungroup() %>% select(!!!grps, !!!groupVars) %>% distinct() %>% group_by(!!!grps) %>% arrange(!!!groupVars) %>% mutate(
    x_int = row_number()
  )
  
  if(max(groupIds$x_int) > 256) stop(paste0("Compression cannot be used on discrete data with more than 256 levels"))
  
  groupIds = groupIds %>% group_by(!!!grps) %>% mutate(
    x_raw = as.raw(x_int-1),
    C_x = n() # the number of non zero classes
  )
  
  tmp = df %>% ungroup() %>% group_by(!!!grps) %>% mutate(
      N = n()
    ) %>% left_join(groupIds, by = groupsJoin)
  
  tmp2 = tmp %>% group_by(!!!grps, C_x, N) %>% group_modify(function(d,g,...) {
    C0 = length(memCompress(as.raw(rep(0,g$N))))
    C1 = C0+as.double(g$N)/log(g$C_x) #max(sapply(c(1:10), function(i) length(memCompress(as.raw(sample.int(g$C_x,size=g$N,replace=TRUE)-1)))))
    C = length(memCompress(as.vector(d$x_raw)))
    if (C > C1) C=C1 # prevent entropy exceeding theoretical maximum
    H = as.double(C-C0)/(C1-C0) * log(g$C_x) # original paper includes a degrees of freedom parameter here. with this setup this can only be one...?
    # TODO: this C_x term is responsible I think for divergence - it is supposed to be max possible entropy but this changes depending on the number of bins populated.
    # browser()
    return(tibble(
      I = H,
      I_sd = NA,
      method = "Compression"
    ))
  })
  
  return(tmp2 %>% ungroup() %>% select(-C_x))
}


# 
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-461
# 
# $$I = H_x + H_y - H_xy$$
# 
# Where I is information
# 
# x is a member of X and |X| is number of possible values of x with no non empty bins (classes_X):
# y is a member of Y and |Y| is number of possible values of y with no non empty bins (classes_Y):
# xy is a member of XY and |XY| is number of possible values of combination of x and y with no non empty bins (classes_XY):
# N is the number of samples
# 
# $$H_x__mm = H_x__emp + (|Z_x| - 1)/2N$$
# 
# $$I__mm = I__emp + (|X| - 1)/2N + (|Y| - 1)/2N - (|XY| - 1)/2N$$
# 
# $$I__mm = I__emp + (|X| + |Y| - |XY| - 1)/2N$$
# 
# This is a Miller-Madow adjustment.
