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
  Rcpp::IntegerMatrix edges; // Matrix of edges with columns parent,child
  Rcpp::NumericVector edgeLength; // Length of edges, indexed by the id# of the child
  int nAlleles; // Number of alleles
  int nTips; // Number of tips on the tree
  int nNode; // Number of nodes on the tree
  int root; // Index of the root node
  
  // Functions
  Rcpp::NumericMatrix postorderMessagePassing(const Rcpp::NumericVector& data, const Rcpp::NumericVector& rateX, 
                                                   const Rcpp::NumericVector& piX);

public:
  // Functions
  phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
            Rcpp::IntegerVector eGroup,Rcpp::List treeInfo);
  double rate(const int child,const Rcpp::NumericVector& siteX);
  arma::vec pi(const Rcpp::NumericVector& siteX);
  arma::mat rateMatrix(const arma::vec & pi,double rate, double branchLength);
  Rcpp::NumericVector siteLL(const Rcpp::NumericMatrix& data, const Rcpp::NumericMatrix& rateX, 
                                        const Rcpp::NumericMatrix& piX);
  Rcpp::DataFrame getRateIndex();
  Rcpp::DataFrame getPiIndex();
  Rcpp::NumericVector getParams();
  void setParams(Rcpp::NumericVector x, Rcpp::IntegerVector index);
};

#endif