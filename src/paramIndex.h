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
  std::vector<std::vector<int>> lookup;
  // Functions
  Rcpp::IntegerVector order(const Rcpp::IntegerVector y, const Rcpp::IntegerVector z);
  
public:
  // Functions
  paramIndex(Rcpp::IntegerVector grp,Rcpp::IntegerVector col, Rcpp::StringVector nm,int start);
  std::vector<int> getIndex(const std::vector<int> grp, const std::vector<int> col,bool expand);
  Rcpp::DataFrame asDF();
  std::vector<std::vector<int>> getLookup();
  std::vector<int> getGroup();
  std::vector<int> getColumn();
  std::vector<std::string> getName();
};

#endif