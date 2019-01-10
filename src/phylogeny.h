#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramIndex.h"

class phylogeny
{
private:
  // Variables
  Rcpp::NumericVector params; // Vector holding all parameters
  paramIndex rateIndex; // Parameter index object for rate parameters
  paramIndex piIndex; // Parameter index object for allelic stationary distribution parameters
  Rcpp::IntegerVector edgeGroup; // Vector that holds group number for each branch, indexed by id# of child
  int nAlleles;

public:
  // Functions
  phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
            Rcpp::IntegerVector eGroup);
  double rate(const int group,const Rcpp::NumericVector& siteX);
  Rcpp::NumericVector pi(const Rcpp::NumericVector& siteX);
  Rcpp::DataFrame getRateIndex();
  Rcpp::DataFrame getPiIndex();
  Rcpp::NumericVector getParams();
};

#endif