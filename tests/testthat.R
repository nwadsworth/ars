library(testthat)
library(ars)

test_check("ars")

test_that("Function generates correctly for a gamma distribution",{
  num_of_samples = 100000
  y <- rgamma(num_of_samples, shape = 10, scale = 3)
  g <- function(x) {dgamma(x, shape = 10, scale =3)}
  x <- ars(g, nsamples = 1000, min = 0, max = Inf , xinit = c(1,50))
  result <- ks.test(x, y)
  expect_gt(result$p.value, expected = 0.1)
})
