#include "paramIndex.h"
#include "phylogeny.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::IntegerVector rGroups, Rcpp::IntegerVector rCol, Rcpp::StringVector rNm,
                     Rcpp::IntegerVector pGroups, Rcpp::IntegerVector pCol, Rcpp::StringVector pNm){
  params = par;
  // rateIndex = paramIndex(rGroups,rCol,rNm,0);
  // piIndex = paramIndex(pGroups,pCol,pNm,rGroups.length());
  // nAlleles = piIndex.getLookup().ncol();
}

double phylogeny::rate(const int group,const Rcpp::NumericVector& siteX){
  Rcpp::NumericVector par = params[rateIndex.getIndex((Rcpp::IntegerVector) group,
                                                      (Rcpp::IntegerVector) Rcpp::seq(0,siteX.length()),true)];
  return(Rcpp::sum(Rcpp::exp(par*siteX)));
}

double phylogeny::pi(const Rcpp::NumericVector& siteX){
  double a = -1.1;
  return(a);
}

RCPP_MODULE(phylogeny) {
  class_<phylogeny>( "phylogeny" )
  .constructor<Rcpp::NumericVector,Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::StringVector,
  Rcpp::IntegerVector, Rcpp::IntegerVector, Rcpp::StringVector>()
  .method("rate", &phylogeny::rate)
  .method("pi", &phylogeny::pi)
  ;
}