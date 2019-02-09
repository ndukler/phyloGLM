// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "phyloGLM_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// expokit_dgpadm
arma::mat expokit_dgpadm(Rcpp::NumericMatrix& mat, double t, bool transpose);
RcppExport SEXP _phyloGLM_expokit_dgpadm(SEXP matSEXP, SEXP tSEXP, SEXP transposeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    rcpp_result_gen = Rcpp::wrap(expokit_dgpadm(mat, t, transpose));
    return rcpp_result_gen;
END_RCPP
}
// logSumExp
double logSumExp(NumericVector x);
RcppExport SEXP _phyloGLM_logSumExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logSumExp(x));
    return rcpp_result_gen;
END_RCPP
}
// matrixToStlXptr
Rcpp::XPtr<std::vector<std::vector<double>>> matrixToStlXptr(Rcpp::NumericMatrix x);
RcppExport SEXP _phyloGLM_matrixToStlXptr(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixToStlXptr(x));
    return rcpp_result_gen;
END_RCPP
}
// stlMatrixSubset
Rcpp::NumericMatrix stlMatrixSubset(SEXP xpsexp, Rcpp::IntegerVector row, Rcpp::IntegerVector col);
RcppExport SEXP _phyloGLM_stlMatrixSubset(SEXP xpsexpSEXP, SEXP rowSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type xpsexp(xpsexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(stlMatrixSubset(xpsexp, row, col));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_paramIndex();
RcppExport SEXP _rcpp_module_boot_phylogeny();

static const R_CallMethodDef CallEntries[] = {
    {"_phyloGLM_expokit_dgpadm", (DL_FUNC) &_phyloGLM_expokit_dgpadm, 3},
    {"_phyloGLM_logSumExp", (DL_FUNC) &_phyloGLM_logSumExp, 1},
    {"_phyloGLM_matrixToStlXptr", (DL_FUNC) &_phyloGLM_matrixToStlXptr, 1},
    {"_phyloGLM_stlMatrixSubset", (DL_FUNC) &_phyloGLM_stlMatrixSubset, 3},
    {"_rcpp_module_boot_paramIndex", (DL_FUNC) &_rcpp_module_boot_paramIndex, 0},
    {"_rcpp_module_boot_phylogeny", (DL_FUNC) &_rcpp_module_boot_phylogeny, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_phyloGLM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
