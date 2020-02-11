
context("Evaluate classifier")

testData = tibble(
  TP = c(100), # pred pos, obs pos
  FP = c(200), # pred pos, obs neg
  FN = c(150), # pred neg, obs pos
  TN = c(550) # pred neg, obs neg
)

testData2 = tibble(
  X1Y1 = c(100), # pred pos, obs pos
  X1 = c(100+150), # obs pos
  Y1 = c(100+200), # pred pos
  total = c(100+200+150+550) # total
)

# from http://onlineconfusionmatrix.com/
expected = list(
  Sensitivity=0.4000,#	TPR = TP / (TP + FN)
  Specificity=0.7333,#	SPC = TN / (FP + TN)
  Precision=0.3333,#	PPV = TP / (TP + FP)
  NPV=0.7857,#	NPV = TN / (TN + FN)
  FPR=0.2667,#	FPR = FP / (FP + TN)
  FDR=0.6667,#	FDR = FP / (FP + TP)
  FNR=0.6000,#	FNR = FN / (FN + TP)
  Accuracy=0.6500,#	ACC = (TP + TN) / (P + N)
  F1=0.3636,#	F1 = 2TP / (2TP + FP + FN)
  MCC=0.1260#	TP*TN - FP*FN / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
)

test_that("probabilities equivalent", {
  tmp1 = testData %>% probabilitiesFromConfusionMatrix(TP,FP,FN,TN)
  tmp2 = testData2 %>% probabilitiesFromCounts(X1Y1,X1,Y1,total)
  expect_equal(tmp1$p_x1y1, tmp2$p_x1y1)
  expect_equal(tmp1$p_x1, tmp2$p_x1)
  expect_equal(tmp1$p_y0, tmp2$p_y0)
})

test_that("info stats functioning",{
  tmp1 = testData %>% probabilitiesFromConfusionMatrix(TP,FP,FN,TN) %>% calculateConfusionMatrixStats()
  expect_equal(tmp1$false_pos_rate,expected$FPR,tolerance=0.0001)
  expect_equal(tmp1$sensitivity,expected$Sensitivity,tolerance=0.0001)
  expect_equal(tmp1$specificity,expected$Specificity,tolerance=0.0001)
  expect_equal(tmp1$f1,expected$F1,tolerance=0.0001)
  # browser()
})