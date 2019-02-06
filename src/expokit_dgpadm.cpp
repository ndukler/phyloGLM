#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

extern"C" {
  void wrapalldmexpv_(int* n,int* m,double* t,double* v,double* w,double* tol,
                      double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
                      int* ia,int* ja,double*a ,int* nz,double* res,int* mxstep,int* flag,double* flag2,double* flag3);

  void wrapalldgexpv_(int* n,int* m,double* t,double* v,double* w,double* tol,
                      double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
                      int* ia,int* ja,double* a,int* nz,double* res,int* mxstep,int* flag,double* flag2,double* flag3);

  void myDMEXPV_(int* n,int* m,double* t,double* v,double* w,double* tol,
                 double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
                 int* ia,int* ja,double* a,int* nz,int* mxstep,int* flag );

  void myDGEXPV_(int* n,int* m,double* t,double* v,double* w,double* tol,
                 double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
                 int* ia,int* ja,double* a,int* nz,int* mxstep,int* flag );
  void wrapdgpadm_(int* ideg,int* m,double* t,double H[],int* ldh,
                   double wsp[],int* lwsp,int ipiv[],int* iexph,int *ns,int *iflag );
}

arma::mat expokit_dgpadm(const arma::mat& mat, double t, const bool transpose) {
  
  // Check if t <= 0
  if(t <= 0){
    Rcpp::stop("\nSTOP ERROR: expokit_dgpadm() must be provided with t > 0  \n");
  }
  
  // Check that matrix is square
  if(!mat.is_square()){
    Rcpp::stop("\nSTOP ERROR: non-square matrix provided to expokit_dgpadm  \n");
  }
  
  // Set ...
  int ideg = 6;
  
  // Order (numrows/numcols) of the matrix
  int m = mat.n_rows;
  
  // output matrix
  //double res[m*m];
  
  // Transpose matrix, if required and flatten
  int dim = m*m;
  double H[dim];
  // Transpose the input if !transpose to compensate for c++ row-major vs fortran column-major input
  if (!transpose){
    int iter = 0;
    for(int i=0; i<m; i++){
      for(int j=0; j<m; j++){
        H[iter]=mat(j,i);
        iter++;
      }
    }
  } else {
    int iter = 0;
    for(int i=0; i<m; i++){
      for(int j=0; j<m; j++){
        H[iter]=mat(i,j);
        iter++;
      }
    }
  }
  
  // (ldh,m):(input) argument matrix
  int ldh = m;
  
  // lwsp = length of wsp, the workspace
  // wsp(lwsp):(workspace/output) lwsp .ge. 4*m*m+ideg+1
  int lwsp = 4*m*m+ideg+1;
  double wsp[lwsp];
  
  // ipiv(m)   : (workspace)
  int ipiv[m];
  
  // iexph:(output) number such that wsp(iexph) points to exp(tH)
  // i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
  int iexph = 0;
  
  // ns:(output) number of scaling-squaring used
  int ns = 0;
  
  // iflag:(output) exit flag
  // 0 - no problem
  // <0 - problem
  int iflag = 0;
  
  //Run the function:
  wrapdgpadm_(&ideg, &m, &t, H, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);
  // Put results back in a matrix
  std::vector<double> temp(lwsp);
  for(int i=0;i<lwsp;i++){
    temp[i]=wsp[i];
  }
  arma::mat out_mat(m,m);
  int iter= iexph-1;
  for(int i=0; i<m; i++){
    for(int j=0; j<m; j++){
      out_mat(j,i)=wsp[iter];
      iter++;
    }
  }
  return(out_mat);
}

// [[Rcpp::export]]
arma::mat expokit_dgpadm(Rcpp::NumericMatrix& mat, double t, bool transpose){
  arma::mat A(mat.nrow(),mat.ncol());
  for(unsigned int i=0;i<A.n_rows;i++){
    for(unsigned int j=0;j<A.n_cols;j++){
      A(i,j)=mat(i,j);
    }
  }
  arma::mat out = expokit_dgpadm(A, t, transpose);
  return(out);
}

