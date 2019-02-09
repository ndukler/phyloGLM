#ifndef WRAPPER_H
#define WRAPPER_H

#include <Rcpp.h>
using namespace Rcpp;

Rcpp::XPtr<std::vector<std::vector<double>>> matrixToStlXptr(Rcpp::NumericMatrix x);

#endif