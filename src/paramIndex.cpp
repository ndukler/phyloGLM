#include "paramIndex.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Order the elements of y and z;
// Order by y unless there's a tie, then order by z.
Rcpp::IntegerVector paramIndex::order(const Rcpp::IntegerVector y, const Rcpp::IntegerVector z) {
  Rcpp::IntegerVector o = Rcpp::seq(0,y.size()-1);
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

paramIndex::paramIndex(Rcpp::IntegerVector grp,Rcpp::IntegerVector col, Rcpp::StringVector nm,int start):
  lookup(Rcpp::max(grp)+1,std::vector<int>(Rcpp::max(col)+1,-1)){
  if(grp.size()!=col.size()){
    Rcpp::stop("The grp and col vectors must be the same length");
  }
  if(grp.size()!=nm.size()){
    Rcpp::stop("The grp and nm vectors must be the same length");
  }
  Rcpp::IntegerVector ord = order(grp,col);
  group=Rcpp::as<std::vector<int>>(grp[ord]);
  column=Rcpp::as<std::vector<int>>(col[ord]);
  name=Rcpp::as<std::vector<std::string>>(nm[ord]);
  idx.resize(grp.size());
  std::iota(idx.begin(),idx.end(),start);
  for(unsigned int i=0; i<group.size();i++){
    lookup[group[i]][column[i]]=i;
  }
}

std::vector<int> paramIndex::getIndex(int a){
  return(idx);
}

std::vector<int> paramIndex::getIndex(const std::vector<int> grp, const std::vector<int> col,bool expand){
  // Rcpp::Rcout << "Lookup Matrix" << lookup << std::endl;
  // Rcpp::Rcout << "Group" << grp << std::endl;
  // Rcpp::Rcout << "Columns" << col << std::endl;
  std::vector<int> out;
  if(!expand){
    if(grp.size() != col.size()){
      Rcpp::stop("grp.size()!=col.length");
    } 
    else{
      out.resize(grp.size());
      for(unsigned int i=0;i<grp.size();i++){
        int lIdx=lookup[ grp[i] ][ col[i] ];
        if(lIdx==-1){
          Rcpp::stop("Invalid index query.");    
        }
        out[i]=idx[lIdx];
      }
    }
  } 
  else {
    out.resize(grp.size()*col.size());
    int iter=0;
    for(unsigned int i=0;i<grp.size();i++){
      for(unsigned int j=0;j<col.size();j++){
        int lIdx=lookup[grp[i]][col[j]];
        if(lIdx==-1){
          Rcpp::stop("Invalid index query.");    
        }
        out[iter]=idx[lIdx];
        iter++;
      }
    }
  }
  return(out);
}

Rcpp::DataFrame paramIndex::asDF(){
  Rcpp::IntegerVector tempG(group.size());
  Rcpp::IntegerVector tempC(group.size());
  Rcpp::CharacterVector tempN(group.size());
  Rcpp::IntegerVector tempI(group.size());

  //Rcpp::Rcout << tempG << std::endl;
  for(unsigned int i=0; i< group.size();i++){
    tempG(i)=group[i];
    tempC(i)=column[i];
    tempN(i)=name[i];
    tempI(i)=idx[i];
  }
  return(Rcpp::DataFrame::create(_["group"]= tempG, _["column"]= tempC,_["name"]=tempN, _["idx"]=tempI));
  // return(Rcpp::DataFrame::create(_["group"]=Rcpp::IntegerVector(2)));
}

/*
 * Accessor methods
 */
std::vector<std::vector<int>> paramIndex::getLookup(){
  return(lookup);
}

std::vector<int> paramIndex::getGroup(){
  return(group);
}

std::vector<int> paramIndex::getColumn(){
  return(column);
}

std::vector<std::string> paramIndex::getName(){
  return(name);
}

RCPP_MODULE(paramIndex) {
  // helping the compiler disambiguate things
  std::vector<int> (paramIndex::*getIndex_1)(int) = &paramIndex::getIndex;
  std::vector<int> (paramIndex::*getIndex_2)(std::vector<int>,std::vector<int>,bool) = &paramIndex::getIndex;
  
  class_<paramIndex>( "paramIndex" )
  .constructor<IntegerVector,IntegerVector,StringVector,int>()
  .method("asDF", &paramIndex::asDF)
  .method( "getIndex" , getIndex_1)  
  .method( "getIndex" , getIndex_2)
  .method("getLookup", &paramIndex::getLookup)
  .method("getGroup", &paramIndex::getGroup)
  .method("getColumn", &paramIndex::getColumn)
  .method("getName", &paramIndex::getName)
  ;
}