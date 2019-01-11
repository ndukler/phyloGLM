#ifndef LOGSUMEXP_H
#define LOGSUMEXP_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double logSumExp(NumericVector x);
double logSumExpArma(arma::vec x);

#endif