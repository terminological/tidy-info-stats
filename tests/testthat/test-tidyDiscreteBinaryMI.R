context("Missing values MI")

describe("information content of missing features", {
  
  testData = missingData()
  data = testData$data  %>% group_by(feature)
  equivData = testData$equivData 
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = ":memory:")
  dataSQL = con %>% copy_to(data)

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
  
  it("produces same result in SQL", {
    tmp = data %>% group_by(feature) %>% calculateDiscreteAbsentValuesMI(vars(outcome), vars(sample))
    tmp2 = dataSQL %>% group_by(feature) %>% calculateDiscreteAbsentValuesMI(vars(outcome), vars(sample)) %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  })
  
  it("produces same result in SQL", {
    tmp = data %>% group_by(feature) %>% calculateDiscretePresentValuesMI(vars(outcome), vars(sample))
    tmp2 = dataSQL %>% group_by(feature) %>% calculateDiscretePresentValuesMI(vars(outcome), vars(sample)) %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  })
  
  #TODO: check when there are multiple values for every a given sample id (same feature) we get same as if there is none.
  # N.b. use the value as a count for new observations
  
  #TODO: check if this works when there is a fixed sample count or a supplied DF (expectedness?)
  
  
})
