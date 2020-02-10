context("DBPlyr consistency")
message("Tests started")
devtools::load_all("~/Git/tidy-info-stats")
library(tidyinfostats)


testData = missingData %>% group_by(feature,subfeature)

describe("expected counts: one outcome per sample",{
  
  expected = testData %>% expectOnePerSample(vars(outcome),vars(sample_id))
  
  it("expects each feature has 10 samples", {
    expect_equal(expected$N, rep(10,10))
  })
  
  it("expects each outcome to have 5 samples", {
    expect_equal(expected$N_x, rep(5,10))
  })

})

describe("observed versus expected",{
  
  obsVsExp = testData %>% observedVersusExpected(vars(outcome), vars(sample_id))
  it("A1 has 10 observed, 5 of which have a bad outcome", {
    expect_equal(obsVsExp %>% filter(subfeature=="A1") %>% pull(N_obs), c(10,10))
    expect_equal(obsVsExp %>% filter(subfeature=="A1" & outcome=="Bad") %>% pull(N_x_obs), 5)
  })
  
})