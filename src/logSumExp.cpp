#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <algorithm>
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

double logSumExp(std::vector<double> x){
    auto a = *std::max_element(x.begin(), x.end());
    double out=0;
    for(unsigned int i=0; i<x.size();i++){
      out+=std::exp(x[i]-a);
    }
    out=a+std::log(out);
    return(out);
}