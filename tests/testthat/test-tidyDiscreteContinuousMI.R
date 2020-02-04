

degenerateData = read.csv("~/Git/tidy-info-stats/tests/degenerateExample.csv")

debugSubset %>% ungroup() %>% calculateDiscreteContinuousMI(vars(sensitivity), okapi_bm25, method = "KWindow")