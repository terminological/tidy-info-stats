

#' calculate mutual information between a discrete value (X) and a discrete value (Y)
#' 
#' @param df - may be grouped, in which case the value is interpreted as different types of continuous variable
#' @param groupXVars - the column(s) of the discrete value (X) quoted by vars(...)
#' @param groupYVars - the column(s) of the discrete value (Y) quoted by vars(...)
#' @param method - the method employed - valid options are "Empirical","MontgomerySmith","Compression","Histogram","Entropy","Grassberger"
#' @param ... - the other parameters are passed onto the implementations
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteDiscreteMI =  function(df, groupXVars, groupYVars, method="Grassberger", ...) {
  switch (method,
          Empirical = calculateDiscreteDiscreteMI_Empirical(df, groupXVars, groupYVars, ...),
          MontgomerySmith = calculateDiscreteDiscreteMI_Entropy(df, groupXVars, groupYVars, entropyMethod="MontgomerySmith", ...),
          Grassberger = calculateDiscreteDiscreteMI_Entropy(df, groupXVars, groupYVars, entropyMethod="Grassberger", ...),
          Compression = calculateDiscreteDiscreteMI_Entropy(df, groupXVars, groupYVars, entropyMethod="Compression", ...),
          Histogram = calculateDiscreteDiscreteMI_Entropy(df, groupXVars, groupYVars, entropyMethod="Histogram", ...),
          Entropy = calculateDiscreteDiscreteMI_Entropy(df, groupXVars, groupYVars, ...),
          {stop(paste0(method," not a valid option"))}
  )
}

#' calculate mutual information between a discrete value (X) and a discrete value (Y)
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupXVars - the column of the discrete value (X) quoted by vars(...)
#' @param groupYVars - the column of the discrete value (Y)
#' @return a dataframe containing the distinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteDiscreteMI_Empirical = function(df, groupXVars, groupYVars, ...) {
  df %>% probabilitiesFromCooccurrence(groupXVars, groupYVars) %>% calculateMultiClassMI()
}

#' calculate mutual information between a discrete value (X) and a discrete value (Y) using estimates of entropy
#' 
#' @param df - may be grouped, in which case the grouping is interpreted as different types of discrete variable
#' @param groupXVars - the column of the discrete value (X) quoted by vars(...)
#' @param groupYVars - the column of the discrete value (Y)
#' @param entropyMethod - the method used to calculate the entropy (see ?tidyinfostats::calculateDiscreteEntropy) - defaults to "Grassberger"
#' @return a dataframe containing the disctinct values of the groups of df, and for each group a mutual information column (I). If df was not grouped this will be a single entry
#' @export
calculateDiscreteDiscreteMI_Entropy = function(df, groupXVars, groupYVars, entropyMethod="Grassberger", ...) {
  grps = df %>% groups()
  joinList = df %>% joinList(defaultJoin = "join")
  
  H_y = df %>% group_by(!!!grps) %>% 
    calculateDiscreteEntropy(groupYVars, method = entropyMethod, ...) %>% 
    rename(I_y = I, I_y_sd = I_sd) %>% mutate(join = 1)
  
  H_y_given_x = df %>% group_by(!!!grps, !!!groupXVars) %>% 
    calculateDiscreteEntropy(groupYVars, method = entropyMethod, ...) %>% 
    rename(N_x = N, I_given_x = I, I_given_x_sd = I_sd) %>% mutate(join = 1)
  
  suppressWarnings({
    tmp2 = H_y_given_x %>% left_join(H_y, by=joinList) %>% mutate(p_x = as.double(N_x)/N) %>% group_by(!!!grps, N, I_y, I_y_sd) %>% summarise(
      na_check = sum(ifelse(is.na(I_given_x),1,0)),
      I_given_x = sum(I_given_x*p_x,na.rm = TRUE), 
			I_given_x_sd = sum(I_given_x_sd*p_x,na.rm = TRUE)
    	) # sometimes max term has no non NA values. In which case this is NAl, which is ok but generates a warning
  })
  
  tmp2 = tmp2 %>% mutate(
    I = ifelse(na_check == 0 & I_given_x < I_y, I_y - I_given_x, NA),
    I_sd = ifelse(na_check == 0 & I_given_x < I_y, I_y_sd + I_given_x_sd, NA),
    method =  paste0("Entropy - ",entropyMethod)
  ) %>% ungroup() %>% select(!!!grps, N, I, I_sd, method) %>% mutate(I = as.double(I), I_sd = as.double(I_sd))
  
  # browser()
  return(tmp2)
}


