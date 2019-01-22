#include "phylogeny.h"
#include "paramIndex.h"
#include "logSumExp.h"
#include "expokit.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup, Rcpp::List treeInfo) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()){
  params = par;
  nAlleles = piIndex.getLookup().nrow()+1; // Add one since base "0" allele level has no parameters
  edgeGroup = eGroup;
  edges = as<IntegerMatrix>(treeInfo[0]);
  edgeLength = as<NumericVector>(treeInfo[1]);
  nTips = as<int>(treeInfo[2]);
  nNode = edges.nrow()+1;
  root=edges(edges.nrow()-1,0); // Get the index of the root node
}

/*
 * Phylogenetic computations
 */ 

// Compute rate for a branch with child at a given site
double phylogeny::rate(const int child,const Rcpp::NumericVector& siteX){
  Rcpp::NumericVector par = params[rateIndex.getIndex( Rcpp::IntegerVector::create(edgeGroup(child)),
                                                      (Rcpp::IntegerVector) Rcpp::seq(0,siteX.size()-1),true)];
  if(siteX.size()!=par.size()){
    Rcpp::stop("Error in rate calculation: feature/parameter vector length mismatch");
  }
  return(std::exp(Rcpp::sum(par*siteX)));
}

// Compute allele stationary distribution for a given site design matrix  
arma::vec phylogeny::pi(const Rcpp::NumericVector& piV){
  arma::vec p(nAlleles);
  // Set as e^0
  p(0)=1;
  // Rcpp::Rcout << "piInitial: " << p << std::endl;
  // Rcpp::Rcout << "piV: " << piV << std::endl;
  // Compute the numerators of the softmax function 
  for(int i=1; i < nAlleles; i++){
     // Rcpp::Rcout << "Calc pi for allele: " << i << std::endl;
     // Get softmax params, need to subtract one from group index since there are no parameters for base
     // level "0" of allele
    Rcpp::NumericVector par = params[piIndex.getIndex( Rcpp::IntegerVector::create(i)-1,
                                                         (Rcpp::IntegerVector) Rcpp::seq(0,piV.size()-1),true)];
   // Rcpp::Rcout << "piParams: " << par << std::endl;
    p(i)=std::exp(-1*Rcpp::sum(par*piV));
   // Rcpp::Rcout << "pi(i): " << p(i) << std::endl;
  }
  // Compute sum of partition function
  double Z = arma::sum(p);
  // Create normalized pi vector
  return(p/Z);
}

arma::mat phylogeny::rateMatrix(const arma::vec& pi,double rate, double branchLength) {
  // initalize temp matrix with all ones
  arma::mat temp(nAlleles,nAlleles,arma::fill::ones);
  // Set diagonal elements to zero
  temp.diag().zeros();
  // constuction that guarentees detailed balance: pi_i*q_ij = pi_j*q_ji
  arma::mat piM = arma::diagmat(pi);
  arma::mat Q=temp*piM; 
  Q.diag()= -arma::sum(Q,1); 
  // Standardize rate matrix and scale by branch length and rate
  double norm = 1/sum(-Q.diag()%pi);
  // Rcpp::Rcout << "QNorm: " << norm << std::endl;
  Q*=norm*rate;
  // Rcpp::Rcout << "QMat: " << Q << std::endl;
  // exponentiate rate matrix and return
  return(expokit_dgpadm(Q,branchLength,FALSE));
}

/*
 * Tree algorithms
 */

// Forward message passing for a single site
NumericMatrix phylogeny::postorderMessagePassing(const NumericVector& data, const NumericVector& rateV, 
                                                 const NumericVector& piV) {
  // Initialize the message table
  NumericMatrix poTab(nNode,nAlleles);
  // Compute the stationary distribution
  arma::vec sitePi = pi(piV);
  // Rcpp::Rcout << "pi: " << sitePi << std::endl;
  // Initialize values for the tips of the tree
  for(int n=0;n<nTips;n++){
    for(int a=0;a<nAlleles;a++){
      poTab(n,a)=data((n*nAlleles)+a);
    }
  }
  // Rcpp::Rcout << "poTab: " << poTab << std::endl;
  // Now compute the probability for the interior nodes
  for(int n=0;n<edges.nrow();n++){
    // Rcpp::Rcout << "Calculating edge: " << n << std::endl;
    int parentInd=edges(n,0);
    int childInd=edges(n,1);
    //Compute the rate for edge n
    double r = rate(childInd,rateV);
    // Rcpp::Rcout << "Rate: " << r << std::endl;
    // Compute the log rate matrix for edge n 
    arma::mat logTMat = arma::log(rateMatrix(sitePi,r,edgeLength(childInd)));
    // Rcpp::Rcout << "LogTMat: " << logTMat << std::endl;
    // iterate over all parental alleles
    for(int a=0;a<nAlleles;a++){
      // Iterate over child alleles
      Rcpp::NumericVector paths(nAlleles);
      for(int c = 0; c < nAlleles; c++){
        paths(c)=poTab(childInd,c) + logTMat(a,c);
      }
      poTab(parentInd,a) = poTab(parentInd,a) + logSumExp(paths);
    }
  }
  // Rcpp::Rcout << "poTab Final: " << poTab << std::endl;
  return(poTab);
}

/*
 * Log-likelihood functions
 */
Rcpp::NumericVector phylogeny::siteLL(const Rcpp::NumericMatrix& data, const Rcpp::NumericMatrix& rateX, 
                            const Rcpp::NumericMatrix& piX) {
  int sites=data.nrow();
  NumericVector siteLik(sites); // Numeric vector of the for the logProbability of each site
  // Rcpp::Rcout << "sites: " << sites << std::endl;
  // Rcpp::Rcout << "nNode: " << nNode << std::endl;
  // Rcpp::Rcout << "edges: " << edges << std::endl;
  // loop over sites running the post-order message passing algorithm
  for(int i=0;i<sites;i++){
    arma::vec logPi = arma::log(pi(piX(i,_)));
    // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
    Rcpp::NumericVector rootMes = postorderMessagePassing(data(i,_), rateX(i,_), piX(i,_))(root,_);
    Rcpp::NumericVector temp(nAlleles); 
    for(int a = 0; a<nAlleles;a++){
      temp(a) = rootMes(a)+logPi(a);
    }
    siteLik(i)=logSumExp(temp);
  }
  return(siteLik);
}

/*
 * Getter and setter functions
 */

Rcpp::DataFrame phylogeny::getRateIndex(){
  return(rateIndex.asDF());
}

Rcpp::DataFrame phylogeny::getPiIndex(){
  return(piIndex.asDF());
}

Rcpp::NumericVector phylogeny::getParams(){
  return(params);
}

void phylogeny::setParams(Rcpp::NumericVector x, Rcpp::IntegerVector index){
  if(x.size()!=index.size()){
    Rcpp::stop("Error when setting parameter vector: vector length mismatch");
  }
  params[index]=x;
}

RCPP_MODULE(phylogeny) {
  class_<phylogeny>( "phylogeny" )
  .constructor<Rcpp::NumericVector, Rcpp::DataFrame, Rcpp::DataFrame, Rcpp::IntegerVector,Rcpp::List>()
  .method("rate", &phylogeny::rate)
  .method("pi", &phylogeny::pi)
  .method("rateMatrix",&phylogeny::rateMatrix)
  .method("getRateIndex", &phylogeny::getRateIndex)
  .method("getPiIndex", &phylogeny::getPiIndex)
  .method("siteLL", &phylogeny::siteLL)
  .method("getParams", &phylogeny::getParams)
  .method("setParams", &phylogeny::setParams)
  ;
}