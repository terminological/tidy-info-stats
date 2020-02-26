# TODO: move this to data_raw & use usethis::use_data



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
#' @import ClassifierResult
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
# 
# ```{r}
# # devtools::load_all("..")
# testData = bloodResultsSimulation(1000)$data
# 
# #### Continuous probability estimation ----
# 
# ggplot(
#   testData %>% group_by(feature,outcome) %>% tidyinfostats::probabilitiesFromContinuous(value, method="SGolay"),
#   aes(x=value,y=p_x, colour=outcome)) + geom_point() + facet_wrap(vars(feature))
# 
# # debug(probabilitiesFromContinuous_SGolay)
# # debug(applySGolayFilter)
# 
# ggplot(
#   testData %>% group_by(feature,outcome) %>% tidyinfostats::probabilitiesFromContinuous(value, method="Kernel"),
#   aes(x=value,y=p_x, colour=outcome)) + geom_point() + facet_wrap(vars(feature))
# 
# ```

missingData = function() {
  # start with a defintion for our test data
  # feature A is present in 80% of outcome 1; 20% of outcome 2 - there is information in missingness
  # feature B is present in 10% of outcome 1; 10% of outcome 2 - there is no information in missingness
  # feature C is present in 40% of outcome 1; 20% of outcome 2 - there is information but less than in A
  # feature D is present in 100% of outcome 1; 100% of outcome 2 - not missing / no information
  missingness = tibble(
    feature = c("A","A","B","B","C","C","D","D"),
    outcome = c(1,2,1,2,1,2,1,2),
    presence = c(0.8,0.2,0.1,0.1,0.4,0.2,1,1)
  )
  
  # outcome 1 seen in 60% of cases outcome 2 in 40%
  expectedness = tibble(
    outcome = c(1,2),
    expected = c(60,40)
  )
  
  # generate a complete data set with a random value and missingness flag
  equivData = expectedness %>% left_join(missingness, by="outcome") %>% group_by(feature,outcome,expected,presence) %>% group_modify(function(d,g,..) {
    return(tibble(
      value = sample.int(4,size = g$expected, replace = TRUE),
      status = c(rep("present",round(g$presence*g$expected)),rep("absent",round((1-g$presence)*g$expected)))
    ))
  }) %>% group_by(feature) %>% arrange(outcome) %>% mutate(sample = c(1:100))
  
  # create test data set with missing values
  data = equivData %>% filter(status != "absent")
  
  return(list(missingness= missingness, expectedness = expectedness, data=data,equivData=equivData))
  
}