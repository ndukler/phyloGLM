#include "paramIndex.h"
#include "phylogeny.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::IntegerVector rGroups, Rcpp::IntegerVector rCol, Rcpp::StringVector rNm,
                     Rcpp::IntegerVector pGroups, Rcpp::IntegerVector pCol, Rcpp::StringVector pNm) : 
                     rateIndex(rGroups,rCol,rNm,0),piIndex(pGroups,pCol,pNm,rGroups.size()){
  params = par;
  nAlleles = piIndex.getLookup().ncol();
}

double phylogeny::rate(const int group,const Rcpp::NumericVector& siteX){
  Rcpp::NumericVector par = params[rateIndex.getIndex( Rcpp::IntegerVector::create(group),
                                                      (Rcpp::IntegerVector) Rcpp::seq(0,siteX.size()-1),true)];
  return(std::exp(Rcpp::sum(par*siteX)));
}

double phylogeny::pi(const Rcpp::NumericVector& siteX){
  double a = -1.1;
  return(a);
}

Rcpp::DataFrame phylogeny::getRateIndex(){
  return(rateIndex.asDF());
}

Rcpp::DataFrame phylogeny::getPiIndex(){
  return(piIndex.asDF());
}

RCPP_MODULE(phylogeny) {
  class_<phylogeny>( "phylogeny" )
  .constructor<Rcpp::NumericVector,Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::StringVector,
  Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::StringVector>()
  .method("rate", &phylogeny::rate)
  .method("pi", &phylogeny::pi)
  .method("getRateIndex", &phylogeny::getRateIndex)
  .method("getPiIndex", &phylogeny::getPiIndex)
  ;
}