#ifndef PARAMINDEX_H
#define PARAMINDEX_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class paramIndex
{
private:
  // Variables
  std::vector<int> group;
  std::vector<int> column;
  std::vector<std::string> name;
  std::vector<int> idx;
  std::vector<int> lookup;
  // Functions
  Rcpp::IntegerVector order(std::vector<int> & y, std::vector<int> & z);
  
public:
  // Functions
  paramIndex(Rcpp::IntegerVector grp,Rcpp::IntegerVector col, Rcpp::StringVector nm,int start);
  Rcpp::IntegerVector getIndex(Rcpp::IntegerVector grp, Rcpp::IntegerVector col,bool expand);
  Rcpp::DataFrame asDF();
  Rcpp::IntegerMatrix getLookup();
};

#endif