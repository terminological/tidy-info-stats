library(tidyinfostats)
context("Continuous entropy")
nd = NormalDistribution$new(mean=3,sd=1)

describe("entropy calculation close to theoretical", {
  
  set.seed(101)
  data = nd$sample(10000)
  
  it("works for quantile method",{
    tmp = data %>% calculateContinuousEntropy(x, method = "Quantile") %>% pull(I)
    expect_equal(nd$theoreticalEntropy(), tmp, tolerance=0.05)
  })
  
  it("works for PDF method using Kernels",{
    tmp = data %>% calculateContinuousEntropy(x, method = "PDF", probabilityMethod="Kernel") %>% pull(I)
    expect_equal(nd$theoreticalEntropy(), tmp, tolerance=0.01)
  })
  
  it("works for PDF method using SGolay",{
    tmp = data %>% calculateContinuousEntropy(x, method = "PDF", probabilityMethod="SGolay") %>% pull(I)
    expect_equal(nd$theoreticalEntropy(), tmp, tolerance=0.05)
  })
  
  
  
  it("works for quantile method with small samples",{
    total = 0
    set.seed(101)
    for (i in c(1:100)) {
      data = nd$sample(30)
      tmp = data %>% calculateContinuousEntropy(x, method = "Quantile") %>% pull(I)
      total = total + tmp
    }
    expect_equal(nd$theoreticalEntropy(), total/100, tolerance=0.15)
  })
  
  it("works for PDF method using Kernels with small samples",{
    total = 0
    set.seed(101)
    for (i in c(1:100)) {
      data = nd$sample(30)
      tmp = data %>% calculateContinuousEntropy(x, method = "PDF", probabilityMethod="Kernel") %>% pull(I)
      total = total + tmp
    }
    expect_equal(nd$theoreticalEntropy(), total/100, tolerance=0.2)
  })
  
  it("works for PDF method using SGolay with small samples",{
    total = 0
    set.seed(101)
    for (i in c(1:100)) {
      data = nd$sample(30)
      tmp = data %>% calculateContinuousEntropy(x, method = "PDF", probabilityMethod="SGolay") %>% pull(I)
      total = total + tmp
    }
    expect_equal(nd$theoreticalEntropy(), total/100, tolerance=0.05)
  })
  
  
})
