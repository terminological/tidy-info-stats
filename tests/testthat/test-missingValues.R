context("DBPlyr consistency")
message("Tests started")
devtools::load_all("~/Git/tidy-info-stats")
library(tidyinfostats)

# a sparse set of test data
# there are 10 unique samples
# depending on the feature there are up to 10 observations of each feature
testData = tibble(
  feature = c(rep("A",20),rep("B",10)), # features A or B,
  subfeature = c(rep("A1",10),rep("A2",5),rep("A3",5),rep("B1",5),rep("B2",5)),
  sample_id = c(1:10,4:8,6:10,1,3,5,7,9,2,4,6,8,10), # up to 10 samples
  value = sample.int(5,30,TRUE) # this is a distraction
) %>% mutate(
  outcome = ifelse(sample_id %% 2 == 0, "Good","Bad") # odd sample id are Bad, Even good
)

testData = testData %>% group_by(feature,subfeature)

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