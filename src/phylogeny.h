#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramIndex.h"

class phylogeny
{
private:
  // Variables
  Rcpp::NumericVector params;
  paramIndex rateIndex;
  paramIndex piIndex;
  int nAlleles;

public:
  // Functions
  phylogeny(Rcpp::NumericVector par,Rcpp::IntegerVector rGroups, Rcpp::IntegerVector rCol, Rcpp::StringVector rNm,
            Rcpp::IntegerVector pGroups, Rcpp::IntegerVector pCol, Rcpp::StringVector pNm);
  double rate(const int group,const Rcpp::NumericVector& siteX);
  double pi(const Rcpp::NumericVector& siteX);
};

#endif