context("Evaluate classifier")
message("Tests started")
library(tidyinfostats)

set.seed(101)
ts = tidyIris() %>% collectAsTrainingSet(sample,outcome,feature,value)
s = sample(length(ts$outcome),size=0.8*length(ts$outcome))
xTrain = ts$matrix[s,] %>% sparseMatrixToSparseMCsr()
xTest = ts$matrix[-s,] %>% sparseMatrixToSparseMCsr()
yTrain = ts$outcome[s]
yTest = ts$outcome[-s]

m=LiblineaR::LiblineaR(data=xTrain,target=yTrain,type=0,cost=0.001,wi=c(setosa=2,versicolor=1,virginica=0.5),bias=1,verbose=FALSE)
p=predict(m,xTest,proba=TRUE,decisionValues=TRUE)

describe("for the test classifier", {
  
  cr = ClassifierResult$fromPredictions(p$probabilities,yTest)
  
  expect_iris_levels = function(x) {
    expect_equal(sort(x),sort(as.character(levels(iris$Species))))
  }
  
  describe("multiclass classifier", {
    
    it("produces a list of outcome classes",{
      expect_iris_levels(cr$getClasses())
    })
    
    it("produces a confusion matrix",{skip("not implemented")})
    
    it("produces a list of binary classification comparisons",{
      tmp = cr$getBinaryClassifiers()
      expect_iris_levels(names(tmp))
      expect_true(sapply(tmp,function(x) "BinaryClassifierResult" %in% class(x)) %>% purrr::reduce(`&&`))
    })
    
    it("plots a family of ROC curves",{
      p <- cr$plotRoc()
      expect_true(is.ggplot(p))
      expect_known_output(p,"multiClass.Rdata")
      # browser()
    })
    
  })
  
  describe("for the binary classification",{
  
    cr$getBinaryClassifiers()$setosa
      
    it("plots a roc curve",{
      p <- cr$getBinaryClassifiers()$virginica$plotRoc()
      expect_true(is.ggplot(p))
      expect_known_output(p,"binClass.Rdata")
    })
    
    it("calculates an AUC",{skip("not implemented")})
    it("calculates a bootstrapped AUC confidence interval",{skip("not implemented")})
    it("calculates a mutual information",{skip("not implemented")})  
    
  })
  
  
  
})