#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "paramIndex.h"

class phylogeny
{
private:
  // Variables
  std::vector<double> params; // Vector holding all parameters
  paramIndex rateIndex; // Parameter index object for rate parameters
  paramIndex piIndex; // Parameter index object for allelic stationary distribution parameters
  std::vector<int> edgeGroup; // Vector that holds group number for each branch, indexed by id# of child
  std::vector<std::vector<int>> edges; // Matrix of edges with columns parent,child
  std::vector<double> edgeLength; // Length of edges, indexed by the id# of the child
  int nAlleles; // Number of alleles
  int nTips; // Number of tips on the tree
  int nNode; // Number of nodes on the tree
  int root; // Index of the root node
  
  // Functions
  std::vector<std::vector<double>> postorderMessagePassing(const std::vector<double>& data, 
                                const std::vector<double>& rateV, const std::vector<double>& piV);
  void chunkLL(std::vector<double>& siteLik, const std::vector<std::vector<double>>& data, 
               const std::vector<std::vector<double>>& rateX, const std::vector<std::vector<double>>& piX,
               unsigned int start, unsigned int end);
  void test(std::vector<double>& siteLik,int start, int end);
    

public:
  // Functions
  phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
            Rcpp::IntegerVector eGroup,Rcpp::List treeInfo);
  double rate(const int child,const std::vector<double>& siteX);
  arma::vec pi(const std::vector<double>& piV);
  arma::mat rateMatrix(const arma::vec & pi,const double rate, const double branchLength);
  std::vector<double> siteLL(const Rcpp::NumericMatrix& data, const Rcpp::NumericMatrix& rateX, 
                                        const Rcpp::NumericMatrix& piX, const unsigned int threads);
  Rcpp::DataFrame getRateIndex();
  Rcpp::DataFrame getPiIndex();
  std::vector<double> getParams();
  void setParams(Rcpp::NumericVector x, Rcpp::IntegerVector index);
};

#endif