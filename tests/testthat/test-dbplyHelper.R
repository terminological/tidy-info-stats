context("DBPlyr consistency")
message("Tests started")
library(tidyinfostats)

set.seed(101)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = ":memory:")

testData = bloodResultsSimulation(1000)$data
testDataSQL = con %>% copy_to(testData)

tidyUSA = tidyUSArrests()
tidyUSASQL = con %>% copy_to(tidyUSA)

tidyIris = tidyIris()
tidyIrisSQL = con %>% copy_to(tidyIris)

#### Sparse matrix ----
test_that("collect as sparse matrix matches input", {
  out = tidyUSA %>% collectAsSparseMatrix(sample,feature,value)
  expect_equal(unname(out[,"Murder"]),USArrests$Murder)
})

test_that("outcome vector matches expectation",{
  expect_equal(
    testData %>% collectOutcomeVector(sample,outcome),
    testDataSQL %>% collectOutcomeVector(sample,outcome)
  )
})

test_that("collect as training set replicates iris data set",{
  tmp = tidyIris %>% collectAsTrainingSet(sample, outcome, feature, value)
  expect_equal(
    unname(tmp$matrix[,1]), 
    iris$Petal.Length
  )
  expect_equal(
    tmp$outcome,
    iris$Species
  )
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
  
  it("collect as training set produces same output in SQL",{
    tmp = tidyIris %>% collectAsTrainingSet(sample, outcome, feature, value)
    tmp2 = tidyIrisSQL %>% collectAsTrainingSet(sample, outcome, feature, value)
    expect_equal(tmp$matrix, tmp2$matrix)
    expect_equal(tmp$outcome, tmp2$outcome)
  })
  
  #TODO: move to test discretisation
  it("produces same result for discretisation by value", {
    tmp = testData %>% group_by(feature) %>% tidyinfostats::discretise(value, value_discrete, bins=3, method = "ByValue", noUnicode=TRUE)
    tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::discretise(value, value_discrete, bins=3, method = "ByValue") %>% collect()
    expect_equal(tmp$value_discrete, tmp2$value_discrete)
  })
  
  #TODO: move to test probability
  it("produces same result for PDF estimation - SGolay", {
    tmp = testData %>% group_by(feature,outcome) %>% probabilitiesFromContinuous(value, method="SGolay")
    tmp2 = testDataSQL %>% group_by(feature,outcome) %>% probabilitiesFromContinuous(value, method="SGolay") %>% collect()
    expect_equal(tmp$p_x, tmp2$p_x, tolerance=0.001)
  })
  
})



#### SGolay ----
describe("sgolay calculations match signal", {
  testDf = tibble(x=2*pi*c(0:100)/100) %>% mutate(y=sin(x))
  testDfSQL = con %>% copy_to(testDf)
  
  it("works for local dataframes, differentiation", {
    tmp = testDf %>% applySGolayFilter(y, y_est, k_05=3, p=2, m=1) %>% collect()
    expect_equal(tmp$y_est, signal::sgolayfilt(x = testDf$y, p = 2, n = 7, m = 1, ts=1/101), tolerance=0.00001)
  })

  it("works for remote tables, differentiation", {
    tmp = testDfSQL %>% applySGolayFilter(y, y_est, k_05=3, p=2, m=1) %>% collect()
    expect_equal(tmp$y_est, signal::sgolayfilt(x = testDf$y, p = 2, n = 7, m = 1, ts=1/101), tolerance=0.00001)
  })
  
  it("works for remote tables, smoothing", {
    tmp = testDfSQL %>% applySGolayFilter(y, y_est, k_05=5, p=2, m=0) %>% collect()
    expect_equal(tmp$y_est, signal::sgolayfilt(x = testDf$y, p = 2, n = 11, m = 0, ts=1/101), tolerance=0.00001)
  })
  
})

#### Digamma ----
describe("digamma calculations match base", {
  testDf = tibble(x=c(1:1000))
  testDfSQL = con %>% copy_to(testDf, overwrite=TRUE)
  
  it("works for local dataframes" ,{
    tmp = testDf %>% calculateDigamma(x, dgx)
    expect_equal(tmp$dgx, digamma(tmp$x), tolerance = 0.001)
  })
  
  it("works for remote dataframes" ,{
    tmp = testDfSQL %>% calculateDigamma(x, dgx) %>% collect()
    expect_equal(tmp$dgx, digamma(tmp$x), tolerance = 0.001)
  })
})

