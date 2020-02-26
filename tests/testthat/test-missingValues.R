context("DBPlyr consistency")
message("Tests started")
library(tidyinfostats)


testData = missingData()$data %>% group_by(feature)

describe("expected counts: one outcome per sample",{
  
  expected = testData %>% expectOnePerSample(vars(outcome),vars(sample))
  
  it("expects each feature has 10 samples", {
    expect_equal(expected$N, rep(100,8))
  })
  
  it("expects each outcome to have 60:40 split in samples", {
    expect_equal(expected$N_x, rep(c(60,40),4))
  })
  # browser()
})

describe("observed versus expected",{
  
  obsVsExp = testData %>% observedVersusExpected(vars(outcome), vars(sample))
  it("A1 has 10 observed, 5 of which have a bad outcome", {
    expect_equal(obsVsExp %>% filter(feature=="A") %>% pull(N_obs), c(56,56))
    expect_equal(obsVsExp %>% filter(feature=="A" & outcome==2) %>% pull(N_x_obs), 8)
  })
  
})