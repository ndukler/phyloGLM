// [[Rcpp::depends(RcppArmadillo)]]

#include "phylogeny.h"
#include "paramIndex.h"
#include "logSumExp.h"
#include "expokit.h"
#include <RcppThread.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup, Rcpp::List treeInfo) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()){
  params = par;
  nAlleles = piIndex.getLookup().size()+1; // Add one since base "0" allele level has no parameters
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
  // Initialize range vector 0 .. siteX.size()-1
  std::vector<int> grp(siteX.size());
  std::iota(grp.begin(), grp.end(), 0);
  // Indicies to retrieve from params
  std::vector<int> parInd = rateIndex.getIndex(std::vector<int>(1,edgeGroup(child)),grp,true);
  if(siteX.size() != parInd.size()){
    Rcpp::stop("Error in rate calculation: feature/parameter vector length mismatch");
  }
  double r=0;
  for(unsigned int i=0; i< parInd.size();i++){
    r+= params[parInd[i]]*siteX(i);  
  }
  return(std::exp(r));
}

// Compute allele stationary distribution for a given site design matrix  
arma::vec phylogeny::pi(const Rcpp::NumericVector& piV){
  arma::vec p(nAlleles);
  // Set as e^0
  p(0)=1;
  // Compute the numerators of the softmax function 
  for(int i=1; i < nAlleles; i++){
    // Rcpp::Rcout << "Calc pi for allele: " << i << std::endl;
    // Get softmax params, need to subtract one from group index since there are no parameters for base
    // level "0" of allele
    std::vector<int> grp(piV.size());
    std::iota(grp.begin(), grp.end(), 0);
    // Get parameter indicies
    std::vector<int> parInd = piIndex.getIndex(std::vector<int>(1,i-1),grp,true);
    // Rcpp::Rcout << "piParams: " << par << std::endl;
    int temp=0; // accumulator value for exponenetial
    for(unsigned int j=0;j<parInd.size();j++){
      temp+=params[parInd[j]]*piV(j);
    }
    p(i)=std::exp(-1*temp);
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
  Q*=norm*rate;
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
  // Initialize values for the tips of the tree
  for(int n=0;n<nTips;n++){
    for(int a=0;a<nAlleles;a++){
      poTab(n,a)=data((n*nAlleles)+a);
    }
  }
  // Now compute the probability for the interior nodes
  for(int n=0;n<edges.nrow();n++){
    int parentInd=edges(n,0);
    int childInd=edges(n,1);
    //Compute the rate for edge n
    double r = rate(childInd,rateV);
    // Compute the log rate matrix for edge n 
    arma::mat logTMat = arma::log(rateMatrix(sitePi,r,edgeLength(childInd)));
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
  return(poTab);
}

/*
 * Log-likelihood functions
 */

// Function to be called by threads for parallel execution
void phylogeny::chunkLL(Rcpp::NumericVector& siteLik, const Rcpp::NumericMatrix& data, const Rcpp::NumericMatrix& rateX, 
             const Rcpp::NumericMatrix& piX, int start, int end){
  // for(int i=start;i<end;i++){
  //   arma::vec logPi = arma::log(pi(piX(i,_)));
  //   // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
  //   Rcpp::NumericVector rootMes = postorderMessagePassing(data(i,_), rateX(i,_), piX(i,_))(root,_);
  //   Rcpp::NumericVector temp(nAlleles); 
  //   for(int a = 0; a<nAlleles;a++){
  //     temp(a) = rootMes(a)+logPi(a);
  //   }
  //   siteLik(i)=logSumExp(temp);
  // }
}

void phylogeny::test(std::vector<double>& siteLik, int start, int end){
  for(int i=start;i<end;i++){
    siteLik[i]=i;
  }
}


std::vector<double> phylogeny::siteLL(const Rcpp::NumericMatrix& data, const Rcpp::NumericMatrix& rateX, 
                            const Rcpp::NumericMatrix& piX) {
  int sites=data.nrow();
  int nThreads=2;
  int rows=sites / nThreads; // Number of blocks that parallel loop must be executed over
  int extra = sites % nThreads; // remaining rows for last thread
  int start = 0; // each thread does [start..end)
  int end = sites;
   
  std::vector<double> siteLik(sites); // Numeric vector of the for the logProbability of each site
  // std::vector<std::thread> workers; // vector of worker threads
  // // loop over sites running the post-order message passing algorithm
  // 
  // for(int t=0;t<nThreads;t++){
  //   if (t == nThreads-1){ // last thread does extra rows:
  //     end += extra;
  //   }
  //   workers.push_back(std::thread(&phylogeny::test, this, siteLik, 
  //                                 start, end));
  //   // workers.push_back(std::thread(chunkLL, siteLik, data, rateX, piX, start, end));
  //   start = end;
  //   end = start + rows;
  // }
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

void phylogeny::setParams(const Rcpp::NumericVector x, const Rcpp::IntegerVector index){
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