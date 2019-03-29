testthat::context("alleleData")

set.seed(123)

## Create a test data
species <- c("Atest", "B test", "C.test")
aData <- lapply(as.list(species), function(x) matrix(runif(120), ncol = 3))
names(aData) <- species

badData <- lapply(as.list(species), function(x) matrix(runif(120, max = 100), ncol = 3))
names(badData) <- species

## Create log data
logAData <- lapply(aData, log)
names(logAData) <- species

## Create test trees, both one that will fail, and one that will pass
passTree <- ape::unroot(ape::rtree(n = length(species), tip.label = species))
failTree1 <- ape::unroot(ape::rtree(n = length(species), tip.label = c("A", "B", "C")))
failTree2 <- ape::unroot(ape::rtree(n = length(species) + 1, tip.label = c(species, "D")))

## Test that object is created correctly when not given site info
testthat::test_that("alleleData builds correctly without site info", {
  testthat::expect_s4_class(alleleData(data = aData, tree = passTree), "alleleData")
})

## Build object with log probabilities
testthat::test_that("alleleData builds correctly without site info", {
  testthat::expect_s4_class(alleleData(data = logAData, tree = passTree, logProb = TRUE), "alleleData")
})

## Build object with incorrect probabilities
testthat::test_that("alleleData fails to build with invalid probabilities", {
  testthat::expect_error(alleleData(data = badData, tree = passTree), "alleleData")
})

## Test that object fails when passed mis-matched tree and data
testthat::test_that("alleleData fails to build when incorrect tree tip.labels supplied", {
  testthat::expect_error(alleleData(data = aData, tree = failTree1))
})
testthat::test_that("alleleData fails to build when wrong tree supplied", {
  testthat::expect_error(alleleData(data = aData, tree = failTree2))
})
