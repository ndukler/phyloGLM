#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' updateFixed
//' 
//' This function takes an logical vector, a 0-based vector of indices, and an vector
//' of updated logical values to update the first logical vector without triggering
//' a rebuild of any associated S4 object.
//'
//' @param old A logical vector to be updated
//' @param replacement A logical vecotor of replacement values
//' @param index the 0-based indicies to update in the old vector
//' @export
// [[Rcpp::export(rng = false)]]
void updateFixed(Rcpp::LogicalVector old, Rcpp::LogicalVector replacement, Rcpp::IntegerVector index) {
  if( is_true(any(index < 0)) || is_true(any(index >= old.size()))){
    throw std::range_error("Cannot use indicies < 0 or >= old.size()");
  }
  if(replacement.size() != index.size()){
    Rcpp::stop("Index and replacement vectors must be of the same length");
  }
  old[index] = replacement;
}