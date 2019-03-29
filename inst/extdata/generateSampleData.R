library(phyloGLM)
library(data.table)
library(ape)

outdir=dirname(rstudioapi::getActiveDocumentContext()$path)

## Set seed for consistency
set.seed(123)
## Set some parameters for the simulation
species <- c("A", "B", "C", "D", "E")
nAlleles <- 3
nSites <- 1000


## Create test trees, both one that will fail, and one that will pass
tree <- ape::rtree(n = length(species), tip.label = species)
ape::write.tree(tree,file.path(outdir,"sampleTree.nh"))

## Create a covariate data.frame with labels
siteLabels <- data.table::data.table(
  site=1:nSites,
  cre.class = factor(rep(c("Enhancer", "Promoter"), times = c(3 / 4 * nSites, 1 / 4 * nSites))),
  logXpress = rnorm(n = nSites, mean = 0, sd = 1)
)
data.table::fwrite(siteLabels,file.path(outdir,"siteCovariates.txt"),sep = "\t")

## Set the rate parameter formula
rateFormula <- formula(~ cre.class + logXpress + 0)
## Set parameter values
rateParams <- matrix(c(-1, -4.5, -3.3), nrow = 1)
piParams <- matrix(c(0.3, 0.1, 0, 0, -0.1, 0.5), nrow = 2, ncol = 3, byrow = T)
## Rate bounds
rBounds <- c(10^-2, 10)
## Simulate data
aData <- simulateSites(
  tr = tree, covariateTable = siteLabels, rateFormula = rateFormula, rateParams = rateParams,
  piParams = piParams, rateBounds = rBounds
)
out=data.table::rbindlist(lapply(as.list(species), function(x) {
  data.table(site = 1:nrow(aData$data[[x]]),species = x, allele=aData$data[[x]])
}))
fwrite(out,file.path(outdir,"samplePAF.txt"),sep="\t")
