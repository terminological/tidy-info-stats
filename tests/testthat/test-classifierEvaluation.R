context("Evaluate classifier")
message("Tests started")
library(tidyinfostats)


ts = tidyIris() %>% collectAsTrainingSet(sample,outcome,feature,value)
s = sample(length(ts$outcome),size=0.8*length(ts$outcome))
xTrain = ts$matrix[s,] %>% sparseMatrixToSparseMCsr()
xTest = ts$matrix[-s,] %>% sparseMatrixToSparseMCsr()
yTrain = ts$outcome[s]
yTest = ts$outcome[-s]

m=LiblineaR::LiblineaR(data=xTrain,target=yTrain,type=0,cost=1,bias=1,verbose=FALSE)
p=predict(m,xTest,proba=TRUE,decisionValues=TRUE)

describe("for the test classifier", {
  
  cr = ClassifierResult$new(p$probabilities,yTest)
  
  it("produces a list of outcome classes",{})
  it("produces a confusion matrix",{})
  it("produces a list of binary classification comparisons",{})
  
  
  describe("for the binary classification")
  
  it("plots a roc curve",{})
  it("calculates an AUC",{})
  it("calculates a bootstrapped AUC confidence interval",{})
  it("calculates a mutual information")
  
})