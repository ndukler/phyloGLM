library(phyloGLM)
## Setup dataset for testing
## Create tree
tree <- ape::read.tree(text = "((A,B),C);")
tree <- ape::unroot.phylo(tree)
tree <- ape::reorder.phylo(tree, "postorder")
tree$edge.length <- c(0.25, 0.5, 2)

## Settings
species <- c("A", "B", "C")
siteInfo <- data.frame(A = 0.01, B = 0.2)
rateFormula <- formula(~ A + B)
states <- matrix(c(0, 1, 1), nrow = 1)
colnames(states) <- tree$tip.label

## Compute site probabilities using epiAllele and construct alleleData object
aData <- disCharToProb(states, c(0, 1))
ad <- alleleData(data = aData, tree = tree, siteInfo = siteInfo)

## Construct two edge group edgeTable
et <- getEdgeTable(ad)
et[, edgeGroup := c(0, 1, 2)]

sigmoid <- function(z) {
  1 / (1 + exp(-z))
}

##### Begin Tests ######
## -------------------------------------------------------------------------- ##
## construct rateModel
testthat::context("rateModel object can be constructed")
testthat::expect_s4_class({
  suppressWarnings(rateMod <<- rateModel(data = ad, rateFormula = rateFormula, lineageTable = et))
}, class = "rateModel")

setParams(rateMod, rep(0, 12), 0:11)

## If rate model object was not constructed don't trigger any further tests
if (exists("rateMod")) {
  testthat::context("rateModel object getter/setter functions")
  ## Check that parameter vector is the correct length
  testthat::test_that(
    "Parameter vector is retieved correctly",
    testthat::expect_equal(getParams(model = rateMod), rep(0, 12))
  )
  setParams(rateMod, rep(1, length(getParams(rateMod))), 0:(length(getParams(rateMod)) - 1))
  testthat::test_that(
    "Parameter vector is set correctly",
    testthat::expect_equal({
      setParams(rateMod, c(0.1, 0.2, 0.3), 6:8)
      getParams(rateMod)[7:9]
    }, c(0.1, 0.2, 0.3))
  )
  ## Test that correct rate calculations are performed for each group
  rMax <- getRateBounds(rateMod)[2] ## max rate
  rMin <- getRateBounds(rateMod)[1] ## min rate
  rMix0 <- sigmoid(sum(rateMod@rateDM[1, ]))
  rMix1 <- sigmoid(sum(rateMod@rateDM[1, ]))
  rMix2 <- sigmoid(sum(c(0.1, 0.2, 0.3) * rateMod@rateDM[1, ]))
  tr0 <- (1 - rMix0) * rMin + rMix0 * rMax
  tr1 <- (1 - rMix1) * rMin + rMix1 * rMax
  tr2 <- (1 - rMix2) * rMin + rMix2 * rMax
  testthat::test_that(
    "Branch specific rate calculations",
    testthat::expect_equal(
      {
        r0 <- rateMod@phylogeny$rate(0, rateMod@rateDM[1, ])
        r1 <- rateMod@phylogeny$rate(1, rateMod@rateDM[1, ])
        r2 <- rateMod@phylogeny$rate(2, rateMod@rateDM[1, ])
        c(r0, r1, r2)
      },
      c(tr0, tr1, tr2)
    )
  )

  ## Check that the stationary distribution, pi is computed correctly
  Z <- c(1, exp(-sum(rateMod@piDM[1, ])))
  piProb <- Z / sum(Z)
  testthat::test_that(
    "Allele stationary distribution calculations",
    testthat::expect_equal(as.numeric(rateMod@phylogeny$pi(rateMod@piDM[1, ])), piProb)
  )


  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel object transition matrix calculations")
  ## Compute transition tables by hand and compare to function computed ones
  qBase <- matrix(c(-piProb[2], piProb[2], piProb[1], -piProb[1]), ncol = 2, nrow = 2, byrow = TRUE)
  qBaseNorm <- qBase / (2 * prod(piProb))
  ## Compute rates
  r0 <- rateMod@phylogeny$rate(0, rateMod@rateDM[1, ])
  r1 <- rateMod@phylogeny$rate(1, rateMod@rateDM[1, ])
  r2 <- rateMod@phylogeny$rate(2, rateMod@rateDM[1, ])
  ## Exponentiate rate matrix
  qBaseNormE1 <- log(as.matrix(Matrix::expm(qBaseNorm * getTree(ad)$edge.length[1] * r0)))
  qBaseNormE2 <- log(as.matrix(Matrix::expm(qBaseNorm * getTree(ad)$edge.length[2] * r1)))
  qBaseNormE3 <- log(as.matrix(Matrix::expm(qBaseNorm * getTree(ad)$edge.length[3] * r2)))

  ## Test against phylogeny object rate matrix calculations
  testthat::test_that(
    "Allele stationary distribution calculations",
    testthat::expect_equal(
      log(as.numeric(rateMod@phylogeny$rateMatrix(piProb, r0, getTree(ad)$edge.length[1]))),
      as.numeric(qBaseNormE1)
    )
  )
  testthat::test_that(
    "Allele stationary distribution calculations",
    testthat::expect_equal(
      log(as.numeric(rateMod@phylogeny$rateMatrix(piProb, r1, getTree(ad)$edge.length[2]))),
      as.numeric(qBaseNormE2)
    )
  )
  testthat::test_that(
    "Allele stationary distribution calculations",
    testthat::expect_equal(
      log(as.numeric(rateMod@phylogeny$rateMatrix(piProb, r2, getTree(ad)$edge.length[3]))),
      as.numeric(qBaseNormE3)
    )
  )
  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel message passing algorithms")
  ## Compute post-order messages
  pl1 <- exp(qBaseNormE1) %*% matrix(c(1, 0), ncol = 1)
  pl2 <- exp(qBaseNormE2) %*% matrix(c(0, 1), ncol = 1)
  pl3 <- exp(qBaseNormE3) %*% matrix(c(0, 1), ncol = 1)
  logAlphaAll <- c(log(c(1, 0)), log(c(0, 1)), log(c(0, 1)), log(pl1) + log(pl2) + log(pl3))

  ## Compute the beta table by hand
  beta1 <- t(pl2 * pl3 * piProb) %*% exp(qBaseNormE1)
  beta2 <- t(pl1 * pl3 * piProb) %*% exp(qBaseNormE2)
  beta3 <- t(pl1 * pl2 * piProb) %*% exp(qBaseNormE3)
  logBetaAll <- log(c(beta1, beta2, beta3, piProb))
  ## Compute alpha and beta tables with test function
  abTab <- rateMod@phylogeny$testMsgPassing(ad@data@x, rateMod@rateDM@x, rateMod@piDM@x)
  ## Check pre-order message passing algorithm
  testthat::test_that(
    "Check pre-order message passing table",
    testthat::expect_equal(logBetaAll, unlist(abTab$beta))
  )
  testthat::test_that(
    "Check post-order message passing table",
    testthat::expect_equal(logAlphaAll, unlist(abTab$alpha))
  )

  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel logLikelihood calculation")
  ## Compute log-likelihood
  ## Check that the manually computed log-likelihood is equal to the output of the function
  testthat::test_that(
    "The correct log-likelihood is computed - E3",
    testthat::expect_equal(logSumExp(log(pl1) + log(pl2) + log(pl3) + log(piProb)), siteLL(model = rateMod))
  )

  ## -------------------------------------------------------------------------- ##
  testthat::context("Marginal calculations")
  testthat::test_that(
    "Check single node marginals",
    testthat::expect_equal(
      c(1, 0, 0, 1, 0, 1, ((pl1 * pl2 * pl3) * piProb) / sum(((pl1 * pl2 * pl3) * piProb))),
      exp(unlist(marginal(rateMod)))
    )
  )

  ## -------------------------------------------------------------------------- ##
  testthat::context("Gradient calculations")
  g_num <- as.numeric(numDeriv::jacobian(func = phyloGLM:::scaledLL, method = "simple", 
                                         x = getParams(rateMod), model = rateMod, scale = -1, threads = 1))
  g_ana <- phyloGLM:::phyloGrad(model = rateMod, scale = -1)
  testthat::test_that(
    "Check rate gradient calculations",
    testthat::expect_gt(cor(g_num, g_ana), 0.99)
  )
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Expected transitions calculations")
  ## Compute pairwise marginal between ancestor (D) and all tips
  logPairwisePotential_DA <- matrix(c(logAlphaAll[1]+logBetaAll[7]+qBaseNormE1[1,1], ## P(D = 1, A = 1)
                            logAlphaAll[1]+logBetaAll[8]+qBaseNormE1[2,1], ## P(D = 2, A = 1)
                            logAlphaAll[2]+logBetaAll[7]+qBaseNormE1[1,2],
                            logAlphaAll[2]+logBetaAll[8]+qBaseNormE1[2,2]), byrow = F, nrow = 2)
  logPairwiseProb_DA <- exp(logPairwisePotential_DA-logSumExp(as.numeric(logPairwisePotential_DA)))
  mt <- marginalTransitions(rateMod)
  testthat::test_that("Testing pairwise marginal calculations",
                      testthat::expect_equal(as.numeric(mt[[1]]),as.numeric(logPairwiseProb_DA))
  )
  testthat::test_that("Testing consistency of nodewise and edgewise marginals",
                      testthat::expect_equal(Reduce("+",mt), marginalTransitions(rateMod,aggregate = "node")[[1]])
  )
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Fixing parameter values")
  testthat::test_that("Testing fixing parameters",{
                      setFixed(rateMod,c(TRUE,TRUE),2:3)
                      testthat::expect_equal(rateMod@fixed[3:4],c(TRUE,TRUE))
  })
  
  testthat::test_that("Testing fixing parameters, should get error",{
    testthat::expect_error(setFixed(rateMod,c(TRUE,TRUE),2:4),
                           "Index and replacement vectors must be of the same length")
  })
  
  testthat::test_that("Testing fixing parameters, should get error",{
    testthat::expect_error(setFixed(rateMod,c(TRUE),20),
                           "Cannot use indicies < 0 or >= old.size()")
  })
  
  testthat::test_that("Testing fixing parameters, should get error",{
    testthat::expect_error(setFixed(rateMod,c(TRUE),-1),
                           "Cannot use indicies < 0 or >= old.size()")
  })
  
  testthat::test_that("Test that ll calculations are correct when only passing unfixed params",{
    ll_unfixed = phyloGLM:::scaledLL(model = rateMod)
    phyloGLM:::setFixed(rateMod,rep(FALSE,length(rateMod@fixed)),0:(length(rateMod@fixed)-1))
    phyloGLM:::setFixed(rateMod,rep(TRUE,3),c(2,5,7))
    ll_fixed = phyloGLM:::scaledLL(x = getParams(rateMod)[!rateMod@fixed],model = rateMod, index = which(!rateMod@fixed)-1)
    testthat::expect_equal(ll_unfixed,ll_fixed)
    phyloGLM:::setFixed(rateMod,rep(FALSE,length(rateMod@fixed)),0:(length(rateMod@fixed)-1))
  })
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Test packing and unpacking")
  ## Compute pairwise marginal between ancestor (D) and all tips
  phyloGLM:::setFixed(rateMod,rep(FALSE,length(rateMod@fixed)),0:(length(rateMod@fixed)-1))
  phyloGLM:::setFixed(rateMod,rep(TRUE,3),c(2,5,7))
  
  testthat::test_that("Model can be packed",{
    testthat::expect_s4_class(packed_model <<- pack(rateMod),"packedModel")
  })
  
  testthat::test_that("Model can be unpacked",testthat::expect_s4_class({
    suppressWarnings(unpacked_model <<- unpack(packed_model))
  } ,class = "rateModel"))
  
  testthat::test_that("Packed and unpacked models have the same parameter vectors",
                      testthat::expect_equal(getParams(rateMod),getParams(unpacked_model))
  )
  
  testthat::test_that("Packed and unpacked models have the same fixed vectors",
                      testthat::expect_equal(rateMod@fixed,unpacked_model@fixed)
  )
  
  testthat::test_that("Packed and unpacked models have the same trees",
                      testthat::expect_true(ape::all.equal.phylo(getTree(rateMod),getTree(unpacked_model)))
  )
  
  testthat::test_that("Packed and unpacked models have the same pi data",
                      testthat::expect_true(all.equal(rateMod@piDM[1:rateMod@piDM@nrow, ],
                                                      unpacked_model@piDM[1:unpacked_model@piDM@nrow, ])))
  
  testthat::test_that("Packed and unpacked models have the same rate data",
                      testthat::expect_true(all.equal(rateMod@rateDM[1:rateMod@rateDM@nrow, ],
                                                      unpacked_model@rateDM[1:unpacked_model@rateDM@nrow, ])))
  
  testthat::test_that("Packed and unpacked models produce the same LL",
                      testthat::expect_equal(ll(model=rateMod),ll(model = unpacked_model))
  )
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Test copying")
  
  suppressWarnings(rateMod_original <- rateModel(data = ad, rateFormula = rateFormula, lineageTable = et))
  setParams(rateMod_original, rep(1, 12), 0:11)
  suppressWarnings(rateMod_copy <- copy(rateMod_original))

  testthat::test_that("There are no side effects from modifying the parameter vector of the copy",{
    setParams(rateMod_copy, 3, 2)                  
    testthat::expect_true(getParams(rateMod_copy)[3] != getParams(rateMod_original)[3])
  })
  
  testthat::test_that("There are no side effects from modifying the fixed vector of the copy",{
    updateFixed(rateMod_copy@fixed, TRUE, 2)                  
    testthat::expect_true(rateMod_copy@fixed[3] != rateMod_original@fixed[3])
  })
  
  testthat::test_that("There are no side effects from forcing a parameter value to zero",{
    forceParameterZero(rateMod_copy,5,verbose = F)
    testthat::expect_true(rateMod_copy@fixed[6] != rateMod_original@fixed[6] &&
                            getParams(rateMod_copy)[6] != getParams(rateMod_original)[6])
  })
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Test pairwise computation exclude functionality")
  rateMod_import = unpack(readRDS(system.file("extdata", "rateMod.Rds", package="phyloGLM")))
  edge_pair_trans = marginalTransitions(rateMod_import, aggregate = "edge")
  edge_pair_trans_exclude = marginalTransitions(rateMod_import, aggregate = "edge", excludeNodes = 1:100)
  node_pair_trans = marginalTransitions(rateMod_import, aggregate = "node")
  node_pair_trans_exclude = marginalTransitions(rateMod_import, aggregate = "node", excludeEdges = c(1,2))
  
  testthat::test_that("Excluding specific node values works",{
    testthat::expect_equal(Reduce(f = "+",edge_pair_trans_exclude),
                             Reduce(f = "+",node_pair_trans[-1:-100]))
  })

  testthat::test_that("Excluding specific edge values works",{
    testthat::expect_equal(Reduce(f = "+",edge_pair_trans[which(!getTree(rateMod_import)$edge[,2] %in% c(1,2))]),
                           Reduce(f = "+", node_pair_trans_exclude))
  })
  
}
