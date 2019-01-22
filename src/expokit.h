#ifndef EXPOKIT_H
#define EXPOKIT_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat expokit_dgpadm(arma::mat& mat, double t, bool transpose);
arma::mat expokit_dgpadm(Rcpp::NumericMatrix& mat, double t, bool transpose);

#endif