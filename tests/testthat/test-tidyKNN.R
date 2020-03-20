context("Tidy KNN")
message("Tests started")
library(tidyinfostats)

set.seed(101)

#### KNN ----
test_that("k nearest neighbour is self", {
  out = mtcars %>% mutate(id = rownames(mtcars)) %>% tidyinfostats::findKNN(mtcars %>% mutate(target=rownames(mtcars)),id,target,k = 1,matchVars = vars(mpg,cyl,disp,hp) ) 
  # This will fail because of equal rank ordering 
  expect_equal(out$id,out$target)
})

