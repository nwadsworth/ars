context("Test against true distribution")
library(testthat)
library(ars)

test_that("Function generates correctly for a gamma distribution",{
  num_of_samples = 100000
  y <- rgamma(num_of_samples, shape = 10, scale = 3)
  g <- function(x) {dgamma(x, shape = 10, scale =3)}
  x <- ars(g, nsamples = 10000, min = 0, max = Inf , xinit = c(1,50))
  result <- ks.test(x, y)
  expect_gt(result$p.value, expected = 0.1)
})

###Test check_dist
test_that("Test student t distribution is non-log-concave distribution",{
  g <- function(x) {dt(x,df=10)}
  expect_false(check_dist(g, nsamples = 1000, xinit = c(-10,2)))
})

test_that("test that normal distribution is log-concave",{
g <- function(x) {dnorm(x)}
expect_true(check_dist(g, nsamples = 1000,xinit = c(-1,0.5)))
})

test_that("test that gamma distribution is log-concave",{
g <- function(x) {dgamma(x, shape = 10, scale =3)}
expect_true(check_dist(g, nsamples = 1000, min = 0, max = Inf , xinit = c(1,50)))
})

test_that("test that Cauchy distribution is non-log-concave",{
g <- function(x) {dcauchy(x)}
expect_false(check_dist(g, nsamples = 1000, xinit = c(-2,2,10)))
})

test_that("test that F distribution with df1=10 df2 = 5 is non-log-concave",{
g <- function(x) {df(x,df1=10, df2 = 5)}
expect_false(check_dist(g, nsamples = 1000, min=0, xinit = c(3,10)))
})
##Test against true distribution
test_that("Function generates correctly for normal distribution",{
  num_of_samples = 100000
  y <- rnorm(num_of_samples)
  g <- function(x) {dnorm(x)}
  x <- ars(g, nsamples = 10000, xinit = c(-1,1))
  result <- ks.test(x, y)
  expect_gt(result$p.value, expected = 0.1)
})

test_that("Test that ars would not work for t distribution",{
  g <- function(x) {dt(x,df=10)}
  expect_null(ars(g, nsamples = 1000, xinit = c(-10,2)))
})

test_that("Test that ars would not work for Cauchy distribution",{
  g <- function(x) {dcauchy(x)}
  expect_null(ars(g, nsamples = 1000, xinit = c(-2,2,10)))
})

test_that("Function generates correctly for exponential distribution",{
  num_of_samples = 100000
  y <- rexp(num_of_samples,rate=3)
  g <- function(x) {dexp(x,rate=3)}
  x <- ars(g, nsamples = 10000, min=0, xinit = c(1,3,10))
  result <- ks.test(x, y)
  expect_gt(result$p.value, expected = 0.1)
})

#test for invalid inputs
test_that("Test when min>max",{
  g <- function(x) {dnorm(x)}
  expect_error(ars(g, nsamples = 1000, min=1, max=0,xinit = c(-1,1)))
})

test_that("Test when there is only one number in xinit",{
  g <- function(x) {dnorm(x)}
  expect_error(ars(g, nsamples = 1000, xinit = c(1)))
})
