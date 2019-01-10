#include "phylogeny.h"
#include "paramIndex.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()){
  params = par;
  nAlleles = piIndex.getLookup().ncol();
  edgeGroup=eGroup;
}

double phylogeny::rate(const int group,const Rcpp::NumericVector& siteX){
  Rcpp::NumericVector par = params[rateIndex.getIndex( Rcpp::IntegerVector::create(group),
                                                      (Rcpp::IntegerVector) Rcpp::seq(0,siteX.size()-1),true)];
  return(std::exp(Rcpp::sum(par*siteX)));
}

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

//Accessor Functions

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
  .constructor<Rcpp::NumericVector, Rcpp::DataFrame, Rcpp::DataFrame, Rcpp::IntegerVector>()
  .method("rate", &phylogeny::rate)
  .method("pi", &phylogeny::pi)
  .method("getRateIndex", &phylogeny::getRateIndex)
  .method("getPiIndex", &phylogeny::getPiIndex)
  .method("getParams", &phylogeny::getParams)
  ;
}