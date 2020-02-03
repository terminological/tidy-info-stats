library(tidyinfostats)
context("Continuous entropy")
nd = NormalDistribution$new(mean=3,sd=1)

describe("entropy calculation close to theoretical", {
  
  set.seed(101)
  data = nd$sample(10000)
  cutsDf  = data %>% cutsDfFromVector(seq(0,6,length.out = 100))
  data = data %>% discretise(x, x_discrete, method="Manual", cutsDf, factorise=TRUE)
  
  infTheo = data %>% calculateDiscreteEntropy(vars(x_discrete), method = "InfoTheo") %>% pull(I)
  
  # http://thirdorderscientist.org/homoclinic-orbit/2013/5/8/bridging-discrete-and-differential-entropy
  # essentially the adjustment to convert between continuous and discrete is +log(range/bins) - but in this case the range is infinite unless you ignore the first and last bin
  it("works for discrete method InfoTheo",{
    expect_equal(nd$theoreticalEntropy(), infTheo+log(6/100), tolerance=0.05)
  })
  
  #the method employed - valid options are "MontgomerySmith", "Histogram", "Grassberger", "InfoTheo", "Compression"
  
  # options(warn=1)
  # debug(calculateDiscreteEntropy_MontgomerySmith)
  
  it("works for discrete method, MontgomerySmith",{
    tmp = data %>% calculateDiscreteEntropy(vars(x_discrete), method = "MontgomerySmith") %>% pull(I)
    expect_equal(infTheo, tmp, tolerance=0.01)
  })

  it("works for discrete method, Grassberger",{
    tmp = data %>% calculateDiscreteEntropy(vars(x_discrete), method = "Grassberger") %>% pull(I)
    expect_equal(infTheo, tmp, tolerance=0.01)
  })
  
  it("works for discrete method, Histogram",{
    tmp = data %>% calculateDiscreteEntropy(vars(x_discrete), method = "Histogram") %>% pull(I)
    expect_equal(infTheo, tmp, tolerance=0.01)
  })
  
  it("works for discrete method, Compression",{
    skip("this method is non functional")
    tmp = data %>% calculateDiscreteEntropy(vars(x_discrete), method = "Compression") %>% pull(I)
    expect_equal(infTheo, tmp, tolerance=0.01)
  })
  
})
