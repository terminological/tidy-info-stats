
#' tidy dataframe of the USArrests data
#' 
#' @import dplyr
#' @export
tidyUSArrests = function() {
  USArrests %>% 
    mutate(sample = rownames(USArrests)) %>%
    tidyr::pivot_longer(-sample, names_to = "feature", values_to = "value")
}

#' tidy dataframe of the USArrests data with co-occurence of features
#' 
#' @import dplyr
#' @export
tidyUSArrestsCooccurrence = function() {
  lhs = tidyUSArrests() %>% rename(feature1=feature, value1=value)
  rhs = tidyUSArrests() %>% rename(feature2=feature, value2=value)
  return(lhs %>% inner_join(rhs, by="sample"))
}

#' tidy dataframe of the USArrests data
#' 
#' @import dplyr
#' @export
tidyDiscreteUSArrests = function() {
  infotheo::discretize(USArrests) %>% 
    mutate(sample = rownames(USArrests)) %>%
    tidyr::pivot_longer(-sample, names_to = "feature", values_to = "value")
}

#' tidy dataframe of the USArrests data with co-occurence of features
#' 
#' @import dplyr
#' @export
tidyDiscreteUSArrestsCooccurrence = function() {
  lhs = tidyDiscreteUSArrests() %>% rename(feature1=feature, value1=value)
  rhs = tidyDiscreteUSArrests() %>% rename(feature2=feature, value2=value)
  return(lhs %>% inner_join(rhs, by="sample"))
}

#' tidy dataframe of the Iris data with features & outcomes
#' 
#' @import dplyr
#' @export
tidyIris = function() {
  iris %>% 
    mutate(sample = row_number()) %>% 
    rename(
      Sepal_Length = Sepal.Length,
      Sepal_Width = Sepal.Width,
      Petal_Length = Petal.Length,
      Petal_Width = Petal.Width
      ) %>%
    tidyr::pivot_longer(cols=c(Sepal_Length,Sepal_Width,Petal_Length,Petal_Width), names_to = "feature") %>% rename(outcome = Species)
}

#' tidy dataframe of the simulation of blood test results with known distributions for individual outcomes.
#' 
#' @import dplyr
#' @export
bloodResultsSimulation = function(n,seed = 101) {

  set.seed(seed)
  
  hb = ConditionalDistribution$new()
  hb$withDistribution(LogNormalDistribution$new(mode=12,sd=1.3), "asymptomatic")
  hb$withDistribution(LogNormalDistribution$new(mode=8,sd=1.5), "tired")
  hb$withDistribution(LogNormalDistribution$new(mode=4,sd=2), "unwell")
  
  k = ConditionalDistribution$new()
  k$withDistribution(NormalDistribution$new(mean=1,sd=0.5), "unwell")
  k$withDistribution(NormalDistribution$new(mean=2,sd=1), "asymptomatic")
  k$withDistribution(NormalDistribution$new(mean=8,sd=3), "tired")

  mvd = MultivariableDistribution$new()
  mvd$withConditionalDistribution(hb,"haemoglobin")
  mvd$withConditionalDistribution(k,"serum k")
  
  mvd$withClassWeights(list(
    unwell=0.2,
    asymptomatic=0.6,
    tired=0.3
  ))
  
  return(list(
    theoretical = tibble(
      feature=c("haemoglobin","serum k"),
      I=c(hb$theoreticalMI(),k$theoreticalMI()),
      mean=c(hb$theoreticalMean(),k$theoreticalMean()),
      var=c(hb$theoreticalVariance(),k$theoreticalVariance())
    ),
    data = mvd$sample(n),
    plot = mvd$plot,
    sample = mvd$sample
  ))
  
}


