#include <Rcpp.h>
using namespace Rcpp;

//' matrixToStlPointer
//' 
//' This function converts a matrix to a 2D double stl vector and returns a XPtr to a 
//' pointer to the vector
//'
//' @param x A numeric matrix
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<std::vector<std::vector<double>>> matrixToStlXptr(Rcpp::NumericMatrix x){
  // Initialze stl vector with new, reqired b/c it allow for dynamic scoping
  std::vector<std::vector<double>> * out = new std::vector<std::vector<double>>(x.nrow(),std::vector<double>(x.ncol()));
  // Fill vector
  for(int i=0; i<x.nrow();i++){
    for(int j=0;j<x.ncol();j++){
      (*out)[i][j]=x(i,j);
    }
  }
  // Wrap pointer to vector
  Rcpp::XPtr<std::vector<std::vector<double>>> p(out, true);
  return(p);
}

//' @export
// [[Rcpp::export(rng = false)]]
std::vector<std::vector<double>> testPtrUpdate(SEXP xpsexp){
  // Initialze stl vector
  XPtr<std::vector<std::vector<double>>> p(xpsexp);
  std::vector<std::vector<double>> x = *p;
  x[0][0]=1;
  return(x);
}