context("String length")
message("Tests started")
library(tidyinfostats)

set.seed(101)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = ":memory:")

testData = bloodResultsSimulation(1000)
testDataSQL = con %>% copy_to(testData)

tidyUSA = tidyUSArrests()
tidyUSASQL = con %>% copy_to(tidyUSArrests())

#### Sparse matrix ----
test_that("collect as sparse matrix matches input", {
  out = tidyUSA %>% collectAsSparseMatrix(sample,feature,value)
  expect_equal(unname(out[,"Murder"]),USArrests$Murder)
})

### Group mutate ----
test_that("groupMutate behaves like tidy mutate", {
  tmp = tidyUSA %>% group_by(feature) %>% mutate(
    mean = mean(value, na.rm=TRUE), 
    sd = sd(value, na.rm=TRUE)) %>% arrange(sample,feature)
  tmp2 = tidyUSASQL %>% group_by(feature) %>% groupMutate(
    mean = mean(value, na.rm=TRUE), 
    sd = sd(value, na.rm=TRUE)) %>% arrange(sample,feature) %>% collect()
  expect_equal(tmp$mean, tmp2$mean, tolerance=0.001)
  expect_equal(tmp$sd, tmp2$sd, tolerance=0.001)
})

#### Dbplyr  ----
describe("dbplyr results same as local data frames", {
  
  it("produces same result for collectAsSparseMatrix", {
    tmp = tidyUSA %>% collectAsSparseMatrix(sample,feature,value)
    tmp2 = tidyUSASQL %>% collectAsSparseMatrix(sample,feature,value)
    expect_equal(tmp, tmp2)
  })
  
  it("produces same result for KWindow MI", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "KWindow")
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "KWindow") %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
  it("produces same result for SGolay MI", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "SGolay")
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "SGolay") %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
  it("produces same result for discretisation by value", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::discretise(value, value_discrete, bins=3, method = "ByValue")
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::discretise(value, value_discrete, bins=3, method = "ByValue") %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
  it("produces same result for MI by discretisation by rank - histogram", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByRank", bins=100)
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByRank", bins=100) %>% collect()
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
  it("produces same result for MI by discretisation by value - grassberger", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod = "Grassberger")
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod = "Grassberger")
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
  it("produces same result for PDF estimation - SGolay", {
    tmp = testData %>% group_by(feature,outcome) %>% probabilitiesFromContinuous(value, method="SGolay")
    tmp2 = testDataSQL %>% group_by(feature,outcome) %>% probabilitiesFromContinuous(value, method="SGolay")
    expect_equal(tmp$I, tmp2$I, tolerance=0.001)
    expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  })
  
})



#### Entropy ----
test_that("entropy calculations match infotheo", {
  fail()
})

