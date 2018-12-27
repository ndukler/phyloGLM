#include <Rcpp.h>
using namespace Rcpp;

class paramIndex{
  
private:
  Rcpp::IntegerVector group;
  Rcpp::IntegerVector column;
  Rcpp::StringVector name;
  Rcpp::IntegerVector idx;
  Rcpp::IntegerMatrix lookup;
  
  // Order the elements of y and z;
  // Order by y unless there's a tie, then order by z.
  IntegerVector order(IntegerVector y, IntegerVector z) {
    IntegerVector o = Rcpp::seq(0,y.length()-1);
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
  
public:
  paramIndex(IntegerVector grp,IntegerVector col, StringVector nm,int start){
    if(grp.length()!=col.length()){
      Rcpp::stop("The grp and col vectors must be the same length");
    }
    if(grp.length()!=nm.length()){
      Rcpp::stop("The grp and nm vectors must be the same length");
    }
    IntegerVector ord = order(grp,col);
    group=clone(grp)[ord];
    column=clone(col)[ord];
    name=clone(nm)[ord];
    idx=seq(0,group.length()-1);
    // Create a matrix lookup to give row for each group/column combo
    lookup=Rcpp::IntegerMatrix(max(group),max(column));
    for(int i=0; i<group.length();i++){
      lookup(group(i),column(i))=idx(i);
    }
  }
  
  int getIndex(int group, int column){
    return(idx(lookup(group,column)));
  }
  
  Rcpp::DataFrame asDF(){
    return(DataFrame::create(_["group"]= group, _["column"]= column,_["name"]=name, _["idx"]=idx));
  }  
  
};

RCPP_MODULE(paramIndex) {
  class_<paramIndex>( "paramIndex" )
  .constructor<IntegerVector,IntegerVector,StringVector,int>()
  .method("asDF", &paramIndex::asDF)
  .method("getIndex", &paramIndex::getIndex)
  ;
}