// [[Rcpp::depends(RcppArmadillo)]]

#include "phylogeny.h"
#include "paramIndex.h"
#include "logSumExp.h"
#include "expokit.h"
#include <thread>
#include <RcppArmadillo.h>
#include <mutex>
using namespace Rcpp;

// std::mutex mu;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup, Rcpp::List treeInfo) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()),
                     edges(Rcpp::as<IntegerMatrix>(treeInfo[0]).nrow(), std::vector<int>(2,0)){
  params = Rcpp::as<std::vector<double>>(par);
  nAlleles = piIndex.getLookup().size()+1; // Add one since base "0" allele level has no parameters
  edgeGroup = Rcpp::as<std::vector<int>>(eGroup);
  // Fill edge matrix
  Rcpp::IntegerMatrix tempEMat=Rcpp::as<IntegerMatrix>(treeInfo[0]);
  for(unsigned int i=0;i<edges.size();i++){
    edges[i][0]=tempEMat(i,0);
    edges[i][1]=tempEMat(i,1);
  }
  edgeLength = Rcpp::as<std::vector<double>>(as<Rcpp::NumericVector>(treeInfo[1]));
  nTips = as<int>(treeInfo[2]);
  nNode = edges.size()+1;
  root=edges[edges.size()-1][0]; // Get the index of the root node
}

/*
 * Phylogenetic computations
 */ 

// Compute rate for a branch with child at a given site
double phylogeny::rate(const int child,const std::vector<double>& siteX){
  // Initialize range vector 0 .. siteX.size()-1
  std::vector<int> grp(siteX.size());
  std::iota(grp.begin(), grp.end(), 0);
  // Indicies to retrieve from params
  std::vector<int> parInd = rateIndex.getIndex(std::vector<int>(1,edgeGroup[child]),grp,true);
  if(siteX.size() != parInd.size()){
    Rcpp::stop("Error in rate calculation: feature/parameter vector length mismatch");
  }
  double r=0;
  for(unsigned int i=0; i< parInd.size();i++){
    r+= params[parInd[i]]*siteX[i];  
  }
  return(std::exp(r));
}

// Compute allele stationary distribution for a given site design matrix  
arma::vec phylogeny::pi(const std::vector<double>& piV){
  arma::vec p(nAlleles);
  // Set as e^0
  p(0)=1;
  // Compute the numerators of the softmax function 
  for(int i=1; i < nAlleles; i++){
    // Rcpp::Rcout << "Calc pi for allele: " << i << std::endl;
    // Get softmax params, need to subtract one from group index since there are no parameters for base
    // level "0" of allele
    std::vector<int> col(piV.size());
    std::iota(col.begin(), col.end(), 0);
    // Get parameter indicies
    std::vector<int> parInd = piIndex.getIndex(std::vector<int>(1,i-1),col,true);
    double temp=0; // accumulator value for exponenetial
    for(unsigned int j=0;j<parInd.size();j++){
      //Rcpp::Rcout << j << " " << parInd[j] << " " << params[parInd[j]] << " " << piV[j] << std::endl;
      temp+=params[parInd[j]]*piV[j];
    }
    //Rcpp::Rcout << "temp: " << temp << std::endl;
    p(i)=std::exp(-temp);
    //Rcpp::Rcout << "pi(i): " << p(i) << std::endl;
  }
  // Compute sum of partition function
  double Z = arma::sum(p);
  // Create normalized pi vector
  return(p/Z);
}

arma::mat phylogeny::rateMatrix(const arma::vec& pi,const double rate, const double branchLength) {
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
  // return(Q);
  return(expokit_dgpadm(Q,branchLength,FALSE));
}

/*
 * Tree algorithms
 */

// Forward message passing for a single site
std::vector<std::vector<double>> phylogeny::postorderMessagePassing(const std::vector<double>& data, 
                                                  const std::vector<double>& rateV, const std::vector<double>& piV) {
  // Initialize the message table
  std::vector<std::vector<double>> poTab(nNode,std::vector<double>(nAlleles,0));
  // Compute the stationary distribution
  arma::vec sitePi = pi(piV);
  // Initialize values for the tips of the tree
  for(int n=0;n<nTips;n++){
    for(int a=0;a<nAlleles;a++){
      poTab[n][a]=data[(n*nAlleles)+a];
    }
  }
  // Now compute the probability for the interior nodes
  for(unsigned int n=0;n<edges.size();n++){
    int parentInd=edges[n][0];
    int childInd=edges[n][1];
    //Compute the rate for edge n
    double r = rate(childInd,rateV);
    // Compute the log rate matrix for edge n
    
    arma::mat logTMat = arma::log(rateMatrix(sitePi,r,edgeLength[childInd]));
    // // iterate over all parental alleles
    for(int a=0;a<nAlleles;a++){
      // Iterate over child alleles
      std::vector<double> paths(nAlleles);
      for(int c = 0; c < nAlleles; c++){
        paths[c]=poTab[childInd][c] + logTMat(a,c);
      }
      poTab[parentInd][a] = poTab[parentInd][a] + logSumExp(paths);
    }
  }
  return(poTab);
}

/*
 * Log-likelihood functions
 */

// Function to be called by threads for parallel execution
void phylogeny::chunkLL(std::vector<double>& siteLik, const std::vector<std::vector<double>>& data, 
                        const std::vector<std::vector<double>>& rateX, const std::vector<std::vector<double>>& piX,
                        unsigned int start, unsigned int end){
  for(unsigned int i=start;i<end;i++){
    arma::vec logPi = arma::log(pi(piX[i]));
    // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
    std::vector<double> rootMes = postorderMessagePassing(data[i], rateX[i], piX[i])[root];
    std::vector<double> temp(nAlleles);
    for(int a = 0; a<nAlleles;a++){
      temp[a] = rootMes[a]+logPi[a];
    }
    siteLik[i]=logSumExp(temp);
  }
}

void phylogeny::test(std::vector<double>& siteLik, int start, int end){
  for(int i=start;i<end;i++){
    siteLik[i]=i;
  }
}


std::vector<double> phylogeny::siteLL(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,const unsigned int threads) {
  // Type and dereference external pointers
  XPtr<std::vector<std::vector<double>>> d(dataPtr);
  std::vector<std::vector<double>> data = *d;
  XPtr<std::vector<std::vector<double>>> r(ratePtr);
  std::vector<std::vector<double>> rateX = *r;
  XPtr<std::vector<std::vector<double>>> p(piPtr);
  std::vector<std::vector<double>> piX = *p;
  
  // Math for setting up block size to pass to threads
  unsigned int sites=data.size();
  unsigned int rows=sites / threads; // Number of blocks that parallel loop must be executed over
  unsigned int extra = sites % threads; // remaining rows for last thread
  unsigned int start = 0; // each thread does [start..end)
  unsigned int end = rows;
  
  std::vector<double> siteLik(sites); // Numeric vector of the for the logProbability of each site
  std::vector<std::thread> workers; // vector of worker threads
  // loop over sites running the post-order message passing algorithm
  for(unsigned int t=0;t<threads;t++){
    if (t == threads-1){ // last thread does extra rows:
      end += extra;
    }
    // workers.push_back(std::thread(&phylogeny::test, this, std::ref(siteLik),start, end));
    // chunkLL(std::ref(siteLik), std::ref(dataS), std::ref(rateXS), std::ref(piXS),start, end);
    workers.push_back(std::thread(&phylogeny::chunkLL, this, std::ref(siteLik), std::ref(data),
                                  std::ref(rateX), std::ref(piX),start, end));
    start = end;
    end = start + rows;
  }
  // Wait for all threads to finish
  for (auto& w : workers) {
    w.join();
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

std::vector<double> phylogeny::getParams(){
  return(params);
}

void phylogeny::setParams(const Rcpp::NumericVector x, const Rcpp::IntegerVector index){
  if(x.size()!=index.size()){
    Rcpp::stop("Error when setting parameter vector: vector length mismatch");
  }
  for(unsigned int i=0; i<x.size();i++){
    params[index(i)]=x(i);
  }
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