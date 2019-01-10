#ifndef PARAMINDEX_H
#define PARAMINDEX_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class paramIndex
{
private:
  // Variables
  Rcpp::IntegerVector group;
  Rcpp::IntegerVector column;
  Rcpp::StringVector name;
  Rcpp::IntegerVector idx;
  Rcpp::IntegerMatrix lookup;
  // Functions
  Rcpp::IntegerVector order(Rcpp::IntegerVector y, Rcpp::IntegerVector z);
  
public:
  // Functions
  paramIndex(Rcpp::IntegerVector grp,Rcpp::IntegerVector col, Rcpp::StringVector nm,int start);
  Rcpp::IntegerVector getIndex(Rcpp::IntegerVector grp, Rcpp::IntegerVector col,bool expand);
  Rcpp::DataFrame asDF();
  Rcpp::IntegerMatrix getLookup();
};

#endif