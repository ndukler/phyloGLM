testthat::context("paramIndex")

## Setup variables to create a parameterIndex
x <- c(1, 1, 1, 2, 3, 4)
y <- c(2, 1, 3, 3, 4, 5)
name <- c("a", "b", "c", "d", "e", "f")

testthat::test_that("paramIndex builds correctly", {
  testthat::expect_s4_class(new(phyloGLM:::paramIndex, x, y, name, 3), "Rcpp_paramIndex")
})

## Test that the index is ordered correctly
foo <- new(phyloGLM:::paramIndex, x, y, name, 3)
temp <- data.frame(group = x, column = y, name)
temp <- temp[with(temp, order(group, column)), ]
testthat::test_that("paramIndex sorts correctly", {
  testthat::expect_equal(foo$getName(), as.character(temp$name))
})

## Test that index starts in correct place
testthat::test_that("paramIndex starts at correct value", {
  testthat::expect_equal(foo$asDF()$idx[1], 3)
})

## Test that the index queries correctly
testthat::test_that("paramIndex queries correctly", {
  testthat::expect_equal(foo$getIndex(1, 1, FALSE), 3)
})
testthat::test_that("paramIndex queries correctly", {
  testthat::expect_equal(foo$getIndex(1, 3, FALSE), 5)
})
testthat::test_that("paramIndex queries correctly", {
  testthat::expect_equal(foo$getIndex(4, 5, FALSE), 8)
})
