#include "phylogeny.h"
#include "paramIndex.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup, Rcpp::List treeInfo) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()){
  params = par;
  nAlleles = piIndex.getLookup().ncol();
  edgeGroup = eGroup;
  edges = as<IntegerMatrix>(treeInfo[0]);
  edgeLength = as<NumericVector>(treeInfo[1]);
  nTips = as<int>(treeInfo[2]);
}

/*
 * Phylogenetic computations
 */ 

// Compute rate for a branch with child at a given site
double phylogeny::rate(const int child,const Rcpp::NumericVector& siteX){
  Rcpp::NumericVector par = params[rateIndex.getIndex( Rcpp::IntegerVector::create(edgeGroup(child)),
                                                      (Rcpp::IntegerVector) Rcpp::seq(0,siteX.size()-1),true)];
  return(std::exp(Rcpp::sum(par*siteX)));
}

// Compute allele stationary distribution for a given site design matrix  
Rcpp::NumericVector phylogeny::pi(const Rcpp::NumericVector& siteX){
  NumericVector p(nAlleles);
  // Set as e^0
  p(0)=1;
  // Compute the numerators of the softmax function 
  for(int i=1; i < nAlleles; i++){
    Rcpp::NumericVector par = params[piIndex.getIndex( Rcpp::IntegerVector::create(i),
                                                         (Rcpp::IntegerVector) Rcpp::seq(0,siteX.size()-1),true)];
    p(i)=std::exp(Rcpp::sum(par*siteX));
  }
  // Compute sum of partition function
  double Z = Rcpp::sum(p);
  return(p/Z);
}

/*
 * Tree algorithms
 */

// Forward message passing for a single site
NumericMatrix phylogeny::postorderMessagePassing(const NumericVector& data, const NumericVector& rateX, 
                                                 const NumericVector& piX) {
  NumericMatrix poTab(nNode,nAlleles);
  // Initialize values for the tips of the tree
  for(unsigned int n=0;n<nTips;n++){
    for(unsigned int a=0;a<nAlleles;a++){
      poTab(n,a)=data((n*nAlleles)+a);
    }
  }
  // Now compute the probability for the interior nodes
  for(unsigned int n=0;n<edges.nrow();n++){
    unsigned int parentInd=edges(n,0);
    unsigned int childInd=edges(n,1);
    for(unsigned int a=0;a<nAlleles;a++){ // iterate over all parental alleles
      poTab(parentInd,a) = poTab(parentInd,a) + logSumExp(poTab(childInd,_) + tMat[childInd](a,_));
    }
  }
  return(poTab);
}

/*
 * Accessor Functions
 */

Rcpp::DataFrame phylogeny::getRateIndex(){
  return(rateIndex.asDF());
}

Rcpp::DataFrame phylogeny::getPiIndex(){
  return(piIndex.asDF());
}

Rcpp::NumericVector phylogeny::getParams(){
  return(params);
}

RCPP_MODULE(phylogeny) {
  class_<phylogeny>( "phylogeny" )
  .constructor<Rcpp::NumericVector, Rcpp::DataFrame, Rcpp::DataFrame, Rcpp::IntegerVector,Rcpp::List>()
  .method("rate", &phylogeny::rate)
  .method("pi", &phylogeny::pi)
  .method("getRateIndex", &phylogeny::getRateIndex)
  .method("getPiIndex", &phylogeny::getPiIndex)
  .method("getParams", &phylogeny::getParams)
  ;
}