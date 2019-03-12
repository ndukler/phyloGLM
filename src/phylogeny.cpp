// [[Rcpp::depends(RcppArmadillo)]]

#include "phylogeny.h"
#include "paramIndex.h"
#include "logSumExp.h"
#include "expokit.h"
#include <thread>
#include <mutex>
#include <RcppArmadillo.h>
using namespace Rcpp;

// mutex for use with the marginal transition function
std::mutex marginalTrans;

phylogeny::phylogeny(Rcpp::NumericVector par,Rcpp::DataFrame rDF, Rcpp::DataFrame pDF,
                     Rcpp::IntegerVector eGroup, Rcpp::List treeInfo,Rcpp::List hyper) : 
                     rateIndex(rDF[0],rDF[1],rDF[2],0),piIndex(pDF[0],pDF[1],pDF[2],rDF.nrow()),
                     edges(Rcpp::as<IntegerMatrix>(treeInfo[0]).nrow(), std::vector<int>(2,0)),
                     siblings(Rcpp::as<Rcpp::ListOf<Rcpp::IntegerVector>>(treeInfo[3]).size()){
  params = Rcpp::as<std::vector<double>>(par);
  nAlleles = piIndex.getLookup().size()+1; // Add one since base "0" allele level has no parameters
  edgeGroup = Rcpp::as<std::vector<int>>(eGroup);
  // Fill edge matrix
  Rcpp::IntegerMatrix tempEMat=Rcpp::as<IntegerMatrix>(treeInfo[0]);
  for(unsigned int i=0;i<edges.size();i++){
    edges[i][0]=tempEMat(i,0);
    edges[i][1]=tempEMat(i,1);
  }
  // Fill sibling matrix
  Rcpp::ListOf<Rcpp::IntegerVector> s = Rcpp::as<Rcpp::ListOf<Rcpp::IntegerVector>>(treeInfo[3]);
  for(unsigned int i=0;i<siblings.size();i++){
    //Iterate over siblings
    siblings[i]=Rcpp::as<std::vector<int>>(s[i]);
  }
  edgeLength = Rcpp::as<std::vector<double>>(as<Rcpp::NumericVector>(treeInfo[1]));
  nTips = as<int>(treeInfo[2]);
  nNode = edges.size()+1;
  root=edges[edges.size()-1][0]; // Get the index of the root node
  // Set rate min/max
  Rcpp::NumericVector rateBounds=Rcpp::as<Rcpp::NumericVector>(hyper[0]);
  rMin = rateBounds(0); // min rate
  rMax = rateBounds(1); // max rate
}

/*
 * Phylogenetic computations
 */ 

// Compute rate for a branch with child at a given site, does not scale by branch length
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
  // Guarentee positive r between 0 and 1
  r=1.0/(1.0+std::exp(-r));
  // // Parameterize rate as a mixture distribution
  double rFinal= (1.0-r)*rMin + r*rMax;
  return(rFinal);
}

// Overloaded function for computing rates, can take many sites
std::vector<double> phylogeny::rateV(const int child,std::vector<double> sites,const SEXP ratePtr){
  XPtr<std::vector<std::vector<double>>> r(ratePtr);
  std::vector<std::vector<double>> rateX = *r;
  std::vector <double> out(sites.size()); 
  for(unsigned int i=0;i<sites.size();i++){
    out[i]=rate(child,rateX[sites[i]]);  
  }
  return(out);
}

// Compute allele stationary distribution for a given site design matrix  
arma::vec phylogeny::pi(const std::vector<double>& piV){
  arma::vec p(nAlleles);
  // Set as log(1)
  p(0)=0;
  // Compute the numerators of the softmax function 
  for(int i=1; i < nAlleles; i++){
    // Rcpp::Rcout << "Calc pi for allele: " << i << std::endl;
    // Get softmax params, need to subtract one from group index since there are no parameters for base
    // level "0" of allele
    std::vector<int> col(piV.size());
    std::iota(col.begin(), col.end(), 0);
    // Get parameter indicies
    std::vector<int> parInd = piIndex.getIndex(std::vector<int>(1,i-1),col,true);
    double temp=0; // accumulator value for -z in e^(-z)
    for(unsigned int j=0;j<parInd.size();j++){
      temp-=params[parInd[j]]*piV[j];
    }
    //Rcpp::Rcout << "temp: " << temp << std::endl;
    p(i)=temp;
    //Rcpp::Rcout << "pi(i): " << p(i) << std::endl;
  }
  // // subtract off the max to prevent overflow
  p=p-arma::max(p);
  p=arma::exp(p);
  // Compute sum of partition function
  double Z = arma::sum(p);
  // Create normalized pi vector
  arma::vec pPrime = p/Z;
  // Check if any of the values too small, and if so replace them
  bool renorm = false;
  for(int i=0; i < nAlleles; i++)
  {
    if(pPrime(i)<=exp(-10)){
      renorm=true;
      pPrime(i)=exp(-10);
    }
  }
  // Renormalize vector if necessary
  if(renorm){
    double zPrime=sum(pPrime);
    pPrime=pPrime/zPrime;
  }
  return(pPrime);
}

// Note that the rows are the parent allele and the columns are the child allele
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
 * Message passing algorithms
 */

// Forward message passing for a single site
// Note that poTab[i][j]  i=node index, j= the allele
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

// Backward message passing on a single site
std::vector<std::vector<double>> phylogeny::preorderMessagePassing(const std::vector<std::vector<double>>& alpha,
                                                                   const std::vector<double>& rateV, 
                                                                   const std::vector<double>& piV) {
  std::vector<std::vector<double>> poTab(nNode,std::vector<double>(nAlleles,0));
  // Initialize the root
  arma::vec sitePi = pi(piV);
  for(int a=0;a<nAlleles;a++){
    poTab[root][a]=std::log(sitePi[a]);
  }
  // Pre-compute the transition matricies and store in vector indexed by child 
  // Note: the root index will be empty, that's ok b/c it should never be called
  std::vector<arma::mat> tMat(nNode,arma::mat(nAlleles,nAlleles));
  for(unsigned int n=0;n<edges.size();n++){
    int childInd=edges[n][1];
    double r = rate(childInd,rateV);
    tMat[childInd]=arma::log(rateMatrix(sitePi,r,edgeLength[childInd]));   
  }
  // Now compute the probability for the interior nodes
  for(int n=edges.size()-1;n>=0;n--){
    int parentInd=edges[n][0];
    int childInd=edges[n][1];
    for(int a=0;a<nAlleles;a++){ // iterate over all alleles of focal node
      std::vector<double> msgHolder(nAlleles);
      for(int b=0;b<nAlleles;b++){ // iterate over all states of parental node summing over all siblings
        double parentContrib = poTab[parentInd][b] + tMat[childInd](b,a); // calculate the parent contribution
        double sibContrib = 0; // holds the sibling contribution
        for(unsigned int s=0;s<siblings[childInd].size();s++){ // iterate over siblings
          std::vector<double> sPartial(nAlleles);
          for(int c = 0; c < nAlleles; c++){ // Iterate over sibling alleles
            sPartial[c]=alpha[siblings[childInd][s]][c]+tMat[siblings[childInd][s]](b,c);
          }
          sibContrib+=logSumExp(sPartial);
        }
        msgHolder[b]=parentContrib+sibContrib;
      }
      poTab[childInd][a]=logSumExp(msgHolder);
    }
  }
  return(poTab);
}

/*
 * Marginal calculations
 */

void phylogeny::chunkMarginal(std::vector<std::vector<std::vector<double>>>& marginal, 
                              const std::vector<std::vector<double>>& data,
                              const std::vector<std::vector<double>>& rateX, 
                              const std::vector<std::vector<double>>& piX,
                        unsigned int start, unsigned int end){
  for(unsigned int i=start;i<end;i++){
    arma::vec logPi = arma::log(pi(piX[i]));
    // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
    std::vector<std::vector<double>> alpha = postorderMessagePassing(data[i], rateX[i], piX[i]);
    std::vector<std::vector<double>> beta = preorderMessagePassing(alpha, rateX[i], piX[i]);
    for(int n=0; n<nNode;n++){
      for(int a=0;a<nAlleles;a++){ // iterate over alleles
        // add contributions from above, below and the stationary dist 
        marginal[i][n][a] = beta[n][a]+alpha[n][a]; 
      }
      // Compute log-partition function
      double Z = logSumExp(marginal[i][n]);
      for(int a=0;a<nAlleles;a++){ // iterate over alleles
        marginal[i][n][a] = marginal[i][n][a]-Z; // Normalize to compute the log-marginal 
      }
    }
  }
}

std::vector<std::vector<std::vector<double>>> phylogeny::marginal(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,const unsigned int threads) {
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
  
  std::vector<std::vector<std::vector<double>>> marginal(sites,std::vector<std::vector<double>>(nNode,std::vector<double>(nAlleles))); // Marginal per site distribution
  std::vector<std::thread> workers; // vector of worker threads
  // loop over sites running the post-order message passing algorithm
  for(unsigned int t=0;t<threads;t++){
    if (t == threads-1){ // last thread does extra rows:
      end += extra;
    }
    workers.push_back(std::thread(&phylogeny::chunkMarginal, this, std::ref(marginal), std::ref(data),
                                  std::ref(rateX), std::ref(piX),start, end));
    start = end;
    end = start + rows;
  }
  // Wait for all threads to finish
  for (auto& w : workers) {
    w.join();
  }
  return(marginal);
}

//
// Edgewise marginal transitions
//

void phylogeny::chunkEdgewiseMarginalTransitions(std::vector<std::vector<std::vector<double>>>& expectedTransitions, 
                              const std::vector<std::vector<double>>& data,
                              const std::vector<std::vector<double>>& rateX, 
                              const std::vector<std::vector<double>>& piX,
                              unsigned int start, unsigned int end){
  for(unsigned int i=start;i<end;i++){
    // Computre the stationary distribution for site i
    arma::vec sitePi = pi(piX[i]);
    // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
    std::vector<std::vector<double>> alpha = postorderMessagePassing(data[i], rateX[i], piX[i]);
    std::vector<std::vector<double>> beta = preorderMessagePassing(alpha, rateX[i], piX[i]);
    for(unsigned int e=0; e<edges.size();e++){
      arma::mat transP(nAlleles,nAlleles,arma::fill::zeros); // The probability of each transition
      int parentInd=edges[e][0];
      int childInd=edges[e][1];
      //Compute the rate for edge n
      double r = rate(childInd,rateX[i]);
      // Compute the log rate matrix for edge n
      arma::mat logTMat = arma::log(rateMatrix(sitePi,r,edgeLength[childInd]));
      // Iterate over all parent child combinations
      for(int a=0;a<nAlleles;a++){ // iterate over parent alleles
        for(int b=0;b<nAlleles;b++){ // iterate over child alleles
          transP(a,b) = beta[parentInd][a]+ logTMat(a,b) + alpha[childInd][b];
        }
      }
      // Compute partition function
      double Z = logSumExpArma(arma::vectorise(transP));
      // Accumulate transition probabilities for site
      marginalTrans.lock();
      for(int a=0;a<nAlleles;a++){ // iterate over parent alleles
        for(int b=0;b<nAlleles;b++){ // iterate over child alleles
          expectedTransitions[e][a][b]+=std::exp(transP(a,b)-Z);
        }
      }
      marginalTrans.unlock();
    }
  }
}

std::vector<std::vector<std::vector<double>>> phylogeny::edgewiseMarginalTransitions(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,const unsigned int threads) {
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
  
  std::vector<std::vector<std::vector<double>>> expectedTransitions(edges.size(),std::vector<std::vector<double>>(nAlleles,std::vector<double>(nAlleles))); // Marginal per site distribution
  std::vector<std::thread> workers; // vector of worker threads
  // loop over sites running the post-order message passing algorithm
  for(unsigned int t=0;t<threads;t++){
    if (t == threads-1){ // last thread does extra rows:
      end += extra;
    }
    workers.push_back(std::thread(&phylogeny::chunkEdgewiseMarginalTransitions, this, std::ref(expectedTransitions), std::ref(data),
                                  std::ref(rateX), std::ref(piX),start, end));
    start = end;
    end = start + rows;
  }
  // Wait for all threads to finish
  for (auto& w : workers) {
    w.join();
  }
  return(expectedTransitions);
}


//
// Nodewise marginal transitions
//
void phylogeny::chunkNodewiseMarginalTransitions(std::vector<std::vector<std::vector<double>>>& expectedTransitions, 
                                                 const std::vector<std::vector<double>>& data,
                                                 const std::vector<std::vector<double>>& rateX, 
                                                 const std::vector<std::vector<double>>& piX,
                                                 unsigned int start, unsigned int end){
  for(unsigned int i=start;i<end;i++){
    // Computre the stationary distribution for site i
    arma::vec sitePi = pi(piX[i]);
    // Rcpp::Rcout << "LogPi: " << logPi << std::endl;
    std::vector<std::vector<double>> alpha = postorderMessagePassing(data[i], rateX[i], piX[i]);
    std::vector<std::vector<double>> beta = preorderMessagePassing(alpha, rateX[i], piX[i]);
    for(unsigned int e=0; e<edges.size();e++){
      arma::mat transP(nAlleles,nAlleles,arma::fill::zeros); // The probability of each transition
      int parentInd=edges[e][0];
      int childInd=edges[e][1];
      //Compute the rate for edge n
      double r = rate(childInd,rateX[i]);
      // Compute the log rate matrix for edge n
      arma::mat logTMat = arma::log(rateMatrix(sitePi,r,edgeLength[childInd]));
      // Iterate over all parent child combinations
      for(int a=0;a<nAlleles;a++){ // iterate over parent alleles
        for(int b=0;b<nAlleles;b++){ // iterate over child alleles
          transP(a,b) = beta[parentInd][a]+ logTMat(a,b) + alpha[childInd][b];
        }
      }
      // Compute partition function
      double Z = logSumExpArma(arma::vectorise(transP));
      // Accumulate transition probabilities for site
      for(int a=0;a<nAlleles;a++){ // iterate over parent alleles
        for(int b=0;b<nAlleles;b++){ // iterate over child alleles
          expectedTransitions[i][a][b]+=std::exp(transP(a,b)-Z);
        }
      }
    }
  }
}

std::vector<std::vector<std::vector<double>>> phylogeny::nodewiseMarginalTransitions(SEXP dataPtr, SEXP ratePtr,SEXP piPtr,const unsigned int threads) {
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
  
  std::vector<std::vector<std::vector<double>>> expectedTransitions(sites,std::vector<std::vector<double>>(nAlleles,std::vector<double>(nAlleles))); // Marginal per site distribution
  std::vector<std::thread> workers; // vector of worker threads
  // loop over sites running the post-order message passing algorithm
  for(unsigned int t=0;t<threads;t++){
    if (t == threads-1){ // last thread does extra rows:
      end += extra;
    }
    workers.push_back(std::thread(&phylogeny::chunkNodewiseMarginalTransitions, this, std::ref(expectedTransitions), std::ref(data),
                                  std::ref(rateX), std::ref(piX),start, end));
    start = end;
    end = start + rows;
  }
  // Wait for all threads to finish
  for (auto& w : workers) {
    w.join();
  }
  return(expectedTransitions);
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


std::vector<double> phylogeny::siteLL(const SEXP dataPtr, const SEXP ratePtr,const SEXP piPtr,
                                      const unsigned int threads) {
  // Type and dereference external pointers
  XPtr<std::vector<std::vector<double>>> d(dataPtr);
  std::vector<std::vector<double>> & data = *d;
  XPtr<std::vector<std::vector<double>>> r(ratePtr);
  std::vector<std::vector<double>> & rateX = *r;
  XPtr<std::vector<std::vector<double>>> p(piPtr);
  std::vector<std::vector<double>> & piX = *p;
  
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

double phylogeny::ll(const SEXP dataPtr, const SEXP ratePtr,const SEXP piPtr, double scale,
                                      const unsigned int threads){
  std::vector<double> sLL = siteLL(dataPtr, ratePtr,piPtr,threads);
  double LL=0;
  for (auto& n : sLL)
    LL += n;
  LL=LL*scale;
  return(LL);
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

std::vector<double> phylogeny::getRateBounds(){
  std::vector<double> v = {rMin,rMax};
  return(v);
}

void phylogeny::setRateBounds(double mn, double mx){
  rMin=mn;
  rMax=mx;
}

void phylogeny::setParams(const Rcpp::NumericVector x, const Rcpp::IntegerVector index){
  if(x.size()!=index.size()){
    Rcpp::stop("Error when setting parameter vector: vector length mismatch");
  }
  for(unsigned int i=0; i<x.size();i++){
    params[index(i)]=x(i);
  }
}

/*
 * Some functions to be called by the test suite
 */
Rcpp::ListOf<std::vector<std::vector<double>>> phylogeny::testMsgPassing(SEXP dataPtr, SEXP ratePtr, SEXP piPtr) {
  // Type and dereference external pointers
  XPtr<std::vector<std::vector<double>>> d(dataPtr);
  std::vector<std::vector<double>> data = *d;
  XPtr<std::vector<std::vector<double>>> r(ratePtr);
  std::vector<std::vector<double>> rateX = *r;
  XPtr<std::vector<std::vector<double>>> p(piPtr);
  std::vector<std::vector<double>> piX = *p;
  
  std::vector<std::vector<double>> alpha = postorderMessagePassing(data[0], rateX[0], piX[0]);
  std::vector<std::vector<double>> beta = preorderMessagePassing(alpha, rateX[0], piX[0]);
  Rcpp::ListOf<std::vector<std::vector<double>>> L = List::create( _["alpha"] = alpha , _["beta"] = beta);
  return(L);
}

/*
 * Validator functions for method dispatch in the Rcpp module
 */
// bool get_numericMatrix_valid(SEXP* args, int nargs){
//   if( nargs != 1 ) return false ;
//   if( TYPEOF(args[0]) != INTSXP ) return false ;
//   return ( LENGTH(args[0]) == 1 ) ;
// }

//.method( "rate" , ( NumericVector (phylogeny::*)(int) )(&phylogeny::get) , &get_int_valid )

RCPP_MODULE(phylogeny) {
  class_<phylogeny>( "phylogeny" )
  .constructor<Rcpp::NumericVector, Rcpp::DataFrame, Rcpp::DataFrame, Rcpp::IntegerVector,Rcpp::List,Rcpp::List>()
  .method("rate", &phylogeny::rate )
  .method("rateV", &phylogeny::rateV)
  .method("pi", &phylogeny::pi)
  .method("rateMatrix",&phylogeny::rateMatrix)
  .method("getRateIndex", &phylogeny::getRateIndex)
  .method("getPiIndex", &phylogeny::getPiIndex)
  .method("siteLL", &phylogeny::siteLL)
  .method("ll", &phylogeny::ll)
  .method("marginal", &phylogeny::marginal)
  .method("edgewiseMarginalTransitions", &phylogeny::edgewiseMarginalTransitions)
  .method("nodewiseMarginalTransitions", &phylogeny::nodewiseMarginalTransitions)
  .method("getParams", &phylogeny::getParams)
  .method("setParams", &phylogeny::setParams)
  .method("getRateBounds", &phylogeny::getRateBounds)
  .method("setRateBounds", &phylogeny::setRateBounds)
  .method("testMsgPassing", &phylogeny::testMsgPassing)
  ;
}