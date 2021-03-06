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
  std::vector<std::vector<int>> siblings; // Matrix of edges with columns parent,child
  std::vector<double> edgeLength; // Length of edges, indexed by the id# of the child
  double rMax = 5; // max rate
  double rMin = 0.001; // min rate
  int nAlleles; // Number of alleles
  int nTips; // Number of tips on the tree
  int nNode; // Number of nodes on the tree
  int root; // Index of the root node
  // Special variables for gradient computations
  double eps; // epsilon for gradient evaluation
  int gradEdgeGroup; // edge group for which dL/dz should be evaluated
  int gradAllele; // allele for which dL/dz should be evaluated
  double epsScale; 
  
  // Functions
  std::vector<std::vector<double>> postorderMessagePassing(const std::vector<double>& data, 
                                const std::vector<double>& rateV, const std::vector<double>& piV);
  std::vector<std::vector<double>> preorderMessagePassing(const std::vector<std::vector<double>>& alpha,
                                                                     const std::vector<double>& rateV, 
                                                                     const std::vector<double>& piV);
  void chunkLL(std::vector<double>& siteLik, const std::vector<std::vector<double>>& data, 
               const std::vector<std::vector<double>>& rateX, const std::vector<std::vector<double>>& piX,
               unsigned int start, unsigned int end);
  void test(std::vector<double>& siteLik,int start, int end);
  void chunkMarginal(std::vector<std::vector<std::vector<double>>>& marginal, 
                                const std::vector<std::vector<double>>& data,
                                const std::vector<std::vector<double>>& rateX, 
                                const std::vector<std::vector<double>>& piX,
                                unsigned int start, unsigned int end);
  void chunkEdgewiseMarginalTransitions(std::vector<std::vector<std::vector<double>>>& expectedTransitions, 
                                           const std::vector<std::vector<double>>& data,
                                           const std::vector<std::vector<double>>& rateX, 
                                           const std::vector<std::vector<double>>& piX,
                                           unsigned int start, unsigned int end,
                                           const std::vector<bool> & include);
  void chunkNodewiseMarginalTransitions(std::vector<std::vector<std::vector<double>>>& expectedTransitions, 
                                                   const std::vector<std::vector<double>>& data,
                                                   const std::vector<std::vector<double>>& rateX, 
                                                   const std::vector<std::vector<double>>& piX,
                                                   unsigned int start, unsigned int end,
                                                   const std::vector<bool> & include);
    

public:
  // Functions
  phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
            Rcpp::IntegerVector eGroup,Rcpp::List treeInfo,Rcpp::List hyper);
  double rate(const int child,const std::vector<double>& siteX);
  std::vector<double> rateV(const int child,std::vector<double> sites,const SEXP ratePtr);
  arma::vec pi(const std::vector<double>& piV);
  Rcpp::NumericMatrix piV_Rcpp(std::vector<double> sites,const SEXP piPtr);
  arma::mat rateMatrix(const arma::vec & pi,const double rate, const double branchLength);
  std::vector<double> siteLL(const SEXP dataPtr, const SEXP ratePtr,const SEXP piPtr,
                             const unsigned int threads);
  double ll(const SEXP dataPtr, const SEXP ratePtr,const SEXP piPtr, double scale,
                       const unsigned int threads);
  std::vector<double> grad(const SEXP dataPtr, const SEXP ratePtr,const SEXP piPtr, double scale,
                         const unsigned int threads);
  std::vector<std::vector<std::vector<double>>> marginal(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,
                                                         const unsigned int threads);
  std::vector<std::vector<std::vector<double>>> edgewiseMarginalTransitions(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,
                                                                            Rcpp::IntegerVector & exclude, 
                                                                            const unsigned int threads);
  std::vector<std::vector<std::vector<double>>> nodewiseMarginalTransitions(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,
                                                                            Rcpp::IntegerVector & exclude,
                                                                            const unsigned int threads);
  // Getters
  Rcpp::DataFrame getRateIndex();
  Rcpp::DataFrame getPiIndex();
  std::vector<double> getParams();
  std::vector<double> getRateBounds();
  // Setters
  void setParams(Rcpp::NumericVector x, Rcpp::IntegerVector index);
  void setRateBounds(double mn, double mx);
  // Functions for test suite purposes only
  Rcpp::ListOf<std::vector<std::vector<double>>> testMsgPassing(SEXP dataPtr, SEXP ratePtr, SEXP piPtr);
};

#endif