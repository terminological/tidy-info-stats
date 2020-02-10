context("Missing values MI")

describe("information content of missing features", {
  
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

  # cross check counts of present and missing values in both data sets  
  # data2 = data %>%  
  #   group_by(feature,outcome) %>% 
  #   summarise(observed=n()) %>% 
  #   left_join(expectedness,by="outcome") %>% 
  #   mutate(missing = expected-observed)
  # 
  # equivData2 = equivData %>% 
  #   group_by(feature,outcome,status) %>% 
  #   summarise(observed=n()) %>% 
  #   tidyr::pivot_wider(names_from=status,values_from = observed) %>%
  #   rename(missing = absent, observed = present) %>% mutate(
  #     missing = ifelse(is.na(missing),0,missing),
  #     expected = missing+observed
  #   )
  
  it("computes MI without throwing an error",{
      expect_silent({
        out = data %>% group_by(feature) %>% calculateDiscreteAbsentValuesMI(vars(outcome), vars(sample))
        out2 = data %>% group_by(feature) %>% calculateDiscretePresentValuesMI(vars(outcome), vars(sample))
      })
  })
  
  # data2 = data2 %>% group_by(feature) %>% mutate(
  #   tot_observed = sum(observed),
  #   tot_expected = sum(expected),
  #   tot_missing = sum(missing),
  #   N = tot_expected
  # )
  
  it("gives same answer as equivalent problem",{
    
    test = bind_rows(
      data %>% group_by(feature) %>% calculateDiscreteAbsentValuesMI(vars(outcome), vars(sample)),
      data %>% group_by(feature) %>% calculateDiscretePresentValuesMI(vars(outcome), vars(sample))
    ) %>% group_by(feature) %>% summarise(I = sum(I)) %>% pull(I)
    
    reference = equivData %>% group_by(feature) %>% calculateDiscreteDiscreteMI_Empirical(vars(outcome), vars(status)) %>% pull(I)
    
    reference2 = equivData %>% group_by(feature) %>% calculateDiscreteDiscreteMI(vars(outcome), vars(status), method="Histogram", mm="FALSE") %>% pull(I)
    
    expect_true(isTRUE(all.equal(test, reference, tolerance=0.00001)))
    expect_true(isTRUE(all.equal(test, reference2, tolerance=0.00001)))
    expect_true(isTRUE(all.equal(reference, reference2, tolerance=0.00001)))
    
  })
  
  #TODO: check when there are multiple values for every a given sample id (same feature) we get same as if there is none.
  # N.b. use the value as a count for new observations
  
  #TODO: check if this works when there is a fixed sample count or a supplied DF (expectedness?)
  
  
})
