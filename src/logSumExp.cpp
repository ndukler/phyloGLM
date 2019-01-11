#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath> 
using namespace Rcpp;

//' Vector logSumExp
//' 
//' This function computes the sum(e^x) of a vector x without leaving log space
//'
//' @param x A numeric vector
//' @export
// [[Rcpp::export(rng = false)]]
double logSumExp(NumericVector x){
  double a=max(x);
  double out=a+log(sum(exp(x-a)));
  return(out);
}

//' Vector logSumExpArma
//'
//' This function computes the sum(e^x) of a vector x without leaving log space
//'
//' @param x A numeric vector
double logSumExpArma(arma::vec x){
  double a=arma::max(x);
  double out=a+log(arma::sum(arma::exp(x-a)));
  return(out);
}