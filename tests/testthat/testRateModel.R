library(phyloGLM)
## Setup dataset for testing
## Create tree
tree=ape::read.tree(text = "((A,B),C);")
tree=ape::unroot.phylo(tree)
tree=ape::reorder.phylo(tree,"postorder")
tree$edge.length=c(0.25,0.5,2)

## Settings
species=c("A","B","C")
siteInfo=data.frame(A=0.01,B=0.2)
rateFormula=formula(~A+B)
states=matrix(c(0,1,1),nrow=1)
colnames(states)=tree$tip.label

## Compute site probabilities using epiAllele and construct alleleData object
aData=disCharToProb(states,c(0,1))
ad=alleleData(data=aData,tree=tree,siteInfo = siteInfo)

## Construct two edge group edgeTable
et=getEdgeTable(ad)
et[,edgeGroup:=c(0,1,2)]

##### Begin Tests ######
## -------------------------------------------------------------------------- ##
## construct rateModel 
testthat::context("rateModel object can be constructed")
testthat::expect_s4_class({suppressWarnings(rateMod<-rateModel(data = ad,rateFormula = rateFormula,lineageTable = et))},class = "rateModel")

## If rate model object was not constructed don't trigger any further tests
if(exists('rateMod')){
  testthat::context("rateModel object getter/setter functions")
  ## Check that parameter vector is the correct length
  testthat::test_that("Parameter vector is retieved correctly",
                      testthat::expect_equal(getParams(obj = rateMod),rep(1,12)))
  testthat::test_that("Parameter vector is set correctly",
                      testthat::expect_equal({
                        setParams(rateMod,c(0.1,0.2,0.3),6:8)
                        getParams(rateMod)[7:9]
                      },c(0.1,0.2,0.3)))
  ## Test that correct rate calculations are performed for each group
  testthat::test_that("Branch specific rate calculations",
                      testthat::expect_equal({
                        r0=rateMod@phylogeny$rate(0,rateMod@rateDM[1,])
                        r1=rateMod@phylogeny$rate(1,rateMod@rateDM[1,])
                        r2=rateMod@phylogeny$rate(2,rateMod@rateDM[1,])
                        c(r0,r1,r2)},
                        c(exp(sum(rateMod@rateDM[1,])),exp(sum(rateMod@rateDM[1,])),
                          exp(sum(c(0.1,0.2,0.3)*rateMod@rateDM[1,])))))
  
  ## Check that the stationary distribution, pi is computed correctly
  Z=c(1,exp(-sum(rateMod@piDM[1,])))
  piProb=Z/sum(Z)
  testthat::test_that("Allele stationary distribution calculations",
                      testthat::expect_equal(as.numeric(rateMod@phylogeny$pi(rateMod@piDM[1,])), piProb))


  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel object transition matrix calculations")
  ## Compute transition tables by hand and compare to function computed ones
  qBase=matrix(c(-piProb[2],piProb[2],piProb[1],-piProb[1]),ncol=2,nrow = 2,byrow = TRUE)
  qBaseNorm=qBase/(2*prod(piProb))
  ## Compute rates
  r0=rateMod@phylogeny$rate(0,rateMod@rateDM[1,])
  r1=rateMod@phylogeny$rate(1,rateMod@rateDM[1,])
  r2=rateMod@phylogeny$rate(2,rateMod@rateDM[1,])
  ## Exponentiate rate matrix
  qBaseNormE1=log(as.matrix(Matrix::expm(qBaseNorm*getTree(ad)$edge.length[1]*r0)))
  qBaseNormE2=log(as.matrix(Matrix::expm(qBaseNorm*getTree(ad)$edge.length[2]*r1)))
  qBaseNormE3=log(as.matrix(Matrix::expm(qBaseNorm*getTree(ad)$edge.length[3]*r2)))

  ## Test against phylogeny object rate matrix calculations
  testthat::test_that("Allele stationary distribution calculations",
                      testthat::expect_equal(log(as.numeric(rateMod@phylogeny$rateMatrix(piProb,r0,getTree(ad)$edge.length[1]))),
                                             as.numeric(qBaseNormE1)))
  testthat::test_that("Allele stationary distribution calculations",
                      testthat::expect_equal(log(as.numeric(rateMod@phylogeny$rateMatrix(piProb,r1,getTree(ad)$edge.length[2]))),
                                             as.numeric(qBaseNormE2)))
  testthat::test_that("Allele stationary distribution calculations",
                      testthat::expect_equal(log(as.numeric(rateMod@phylogeny$rateMatrix(piProb,r2,getTree(ad)$edge.length[3]))),
                                             as.numeric(qBaseNormE3)))
  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel message passing algorithms")
  ## Compute post-order messages   
  pl1=exp(qBaseNormE1) %*% matrix(c(1,0),ncol = 1)
  pl2=exp(qBaseNormE2) %*% matrix(c(0,1),ncol=1)  
  pl3=exp(qBaseNormE3) %*% matrix(c(0,1),ncol=1) 
  
  ## Compute the beta table by hand
  beta1=t(pl2*pl3*piProb) %*% exp(qBaseNormE1)
  beta2=t(pl1*pl3*piProb) %*% exp(qBaseNormE2) 
  beta3=t(pl1*pl2*piProb) %*% exp(qBaseNormE3)
  logBetaAll=log(c(beta1,beta2,beta3,piProb))
  ## Compute alpha and beta tables with test function
  abTab=rateMod@phylogeny$testMsgPassing(ad@data@x,rateMod@rateDM@x,rateMod@piDM@x)
  ## Check pre-order message passing algorithm
  testthat::test_that("Check pre-order message passing table",
                      testthat::expect_equal(logBetaAll,unlist(abTab$beta)))
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("rateModel logLikelihood calculation")
  ## Compute log-likelihood
  ## Check that the manually computed log-likelihood is equal to the output of the function
  testthat::test_that("The correct log-likelihood is computed - E3",
                      testthat::expect_equal(logSumExp(log(pl1)+log(pl2)+log(pl3)+log(piProb)),siteLL(obj=rateMod)))
  
  ## -------------------------------------------------------------------------- ##
  testthat::context("Marginal calculations")
  testthat::test_that("Check single node marginals",
                      testthat::expect_equal(c(1,0,0,1,0,1, ((pl1*pl2*pl3)*piProb)/sum(((pl1*pl2*pl3)*piProb))),
                                             exp(unlist(marginal(rateMod)))))

}

