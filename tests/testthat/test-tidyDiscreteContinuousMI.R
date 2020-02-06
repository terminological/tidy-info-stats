context("discrete continuous MI")
message("Tests started")
library(tidyverse)
library(tidyinfostats)


set.seed(101)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = ":memory:")

test = bloodResultsSimulation(1000)
testData = test$data %>% rename(
  discrete = outcome,
  continuous = value
)
testDataSQL = con %>% copy_to(testData)



# debug(calculateDiscreteContinuousMI_KWindow)

# individual test functions
kwin = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "KWindow", ...)
knn = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "KNN", useKWindow=FALSE, ...)
discByRank = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByRank")
info = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod="Entropy", entropyMethod="InfoTheo")
hist = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod="Histogram")
ms = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod="MontgomerySmith")
grass = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod="Entropy", entropyMethod="Grassberger")
comp1 = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Discretise", discretiseMethod="ByValue", binStrategy = linearBySize(8,4,256), mutualInfoMethod="Entropy", entropyMethod="Compression")
comp2 = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Compression")
sgolay = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Entropy", entropyMethod="PDF", probabilityMethod = "SGolay")
quant = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Entropy", entropyMethod="Quantile")
kernel = function(testData, ...) testData %>% group_by(feature) %>% calculateDiscreteContinuousMI(vars(discrete), continuous, method = "Entropy", entropyMethod="PDF", probabilityMethod = "Kernel")

# test gruops
disc = c(discByRank, info, hist, ms, grass, comp1, comp2)
pdfs = c(sgolay,quant,kernel)
nn = c(kwin,knn)
all = c(nn,pdfs,disc)

reliable = c(kwin,knn,info,hist,grass,sgolay,quant)
poor = c(ms,comp1,comp2)
unreliable = all[!(all %in% c(reliable,poor))]

dbplyrSupported = c(kwin,sgolay,quant,hist,grass)

collectResults = function(testData, algos) {
  
  results = NULL
  results = splice(map(algos, function(algo) {
    set.seed(101)
    return(algo(testData) %>% collect())}
  ))
  return(tibble(test=results) %>% unnest(cols=c(test)))
}


describe("approximates theoretical MI for known distribution",{
  
  
})

#### Small sample ----
describe("handles degenerate data",{

  degenerateData = read.csv("~/Git/tidy-info-stats/tests/degenerateExample.csv") %>% rename(
    feature = note_nlp_concept_id,
    discrete = sensitivity,
    continuous = okapi_bm25
  )
  
  it("produces a result larger than zero",{
    result = collectResults(degenerateData,all)
    expect_true((result %>% pull(I)) > 0 || is.na(result %>% pull(I))) %>% reduce(`&&`)
  })
  
})


#### Average sample ----
describe("handles moderate data", {
  
  it("reliable estimates are within 5% of theoretical values", {
    result = collectResults(testData, reliable)
    theo = test$theoretical
    error = result %>% left_join(theo %>% rename(I_theo = I), by="feature") %>% mutate(error = abs((I-I_theo)/I_theo)) %>% pull(error)
    browser()
    expect_true((error < 0.06) %>% reduce(`&&`))
  })
  
  it("unreliable estimates are within 20% of theoretical values", {
    result = collectResults(testData, unreliable)
    theo = test$theoretical
    error = result %>% left_join(theo %>% rename(I_theo = I), by="feature") %>% mutate(error = (I-I_theo)/I_theo) %>% pull(error)
    browser()
    expect_true((error < 0.2) %>% reduce(`&&`))
  })
  
})

#### Dbplyr  ----
describe("dbplyr results same as local data frames", {
  
  result = collectResults(testData ,dbplyrSupported)
  resultSQL = collectResults(testDataSQL ,dbplyrSupported)
  
  it("produces same result for SQL data", {
    tmp = result$I - resultSQL$I
    expect_true((tmp < 0.001) %>% reduce(`&&`))
  })
  
  # it("produces same result for KWindow MI", {
  #   tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "KWindow")
  #   tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "KWindow") %>% collect()
  #   expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  #   expect_equal(tmp$I_sd, tmp2$I_sd, tolerance=0.001)
  # })
  # 
  # it("produces same result for SGolay MI", {
  #   tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "SGolay")
  #   tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "SGolay") %>% collect()
  #   expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  # })
  # 
  # it("produces same result for MI by discretisation by rank - histogram", {
  #   tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByRank", bins=100, noUnicode=TRUE)
  #   tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByRank", bins=100) %>% collect()
  #   expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  # })
  # 
  # it("produces same result for MI by discretisation by value - grassberger", {
  #   tmp = testData %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod = "Grassberger")
  #   tmp2 = testDataSQL %>% group_by(feature) %>% tidyinfostats::calculateDiscreteContinuousMI(vars(outcome), value, method = "Discretise", discretiseMethod="ByValue", mutualInfoMethod = "Grassberger")
  #   expect_equal(tmp$I, tmp2$I, tolerance=0.001)
  #   
  # })
  

})
