#include "paramIndex.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Order the elements of y and z;
// Order by y unless there's a tie, then order by z.
Rcpp::IntegerVector paramIndex::order(Rcpp::IntegerVector y, Rcpp::IntegerVector z) {
  Rcpp::IntegerVector o = Rcpp::seq(0,y.length()-1);
  // Then sort that vector by the values of y and z
  std::sort(o.begin(), o.end(), [&](int i, int j){
    if ( y(i) == y(j) ) {
      return(z(i) < z(j));
    }
    return(y(i) < y(j));
  });
  // Return order
  return(o);
}
  
  
paramIndex::paramIndex(Rcpp::IntegerVector grp,Rcpp::IntegerVector col, Rcpp::StringVector nm,int start){
  if(grp.length()!=col.length()){
    Rcpp::stop("The grp and col vectors must be the same length");
  }
  if(grp.length()!=nm.length()){
    Rcpp::stop("The grp and nm vectors must be the same length");
  }
  Rcpp::IntegerVector ord = order(grp,col);
  group=clone(grp)[ord];
  column=clone(col)[ord];
  name=clone(nm)[ord];
  idx=seq(0,group.length()-1);
  // Create a matrix lookup to give row for each group/column combo
  lookup=Rcpp::IntegerMatrix(max(group)+1,max(column)+1);
  for(int i=0; i<group.length();i++){
    lookup(group(i),column(i))=idx(i);
  }
}

Rcpp::IntegerVector paramIndex::getIndex(Rcpp::IntegerVector grp, Rcpp::IntegerVector col,bool expand){
  Rcpp::IntegerVector out;
  if(!expand){
    if(grp.length() != col.length()){
      stop("grp.length()!=col.length");
    } 
    else{
      out = Rcpp::IntegerVector(grp.length());
      for(int i=0;i<grp.length();i++){
        out(i)=idx(lookup(grp(i),col(i)));
      }
    }
  } 
  else {
    out = Rcpp::IntegerVector(grp.length()*col.length());
    int iter=0;
    for(int i=0;i<grp.length();i++){
      for(int j=0;j<col.length();j++){
        out(iter)=idx(lookup(grp(i),col(j)));
        iter=iter+1;
      }
    }
  }
  return(out);
}

Rcpp::DataFrame paramIndex::asDF(){
  return(Rcpp::DataFrame::create(_["group"]= group, _["column"]= column,_["name"]=name, _["idx"]=idx));
}

Rcpp::IntegerMatrix paramIndex::getLookup(){
  return(lookup);
}

RCPP_MODULE(paramIndex) {
  class_<paramIndex>( "paramIndex" )
  .constructor<IntegerVector,IntegerVector,StringVector,int>()
  .method("asDF", &paramIndex::asDF)
  .method("getIndex", &paramIndex::getIndex)
  .method("getLookup", &paramIndex::getLookup)
  ;
}