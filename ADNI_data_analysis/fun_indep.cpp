#include <RcppEigen.h>
#include <Rcpp.h>
//#include <RcppArmadillo.h>

#include <Eigen/Core>
#include <Eigen/Dense>

//#include <armadillo>

using namespace Rcpp;
using namespace Eigen;
using Eigen::VectorXi;

//using namespace arma;
// This is matrix side calculation

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
// each column of matIY[,i] is weight centered at X_i
Eigen::MatrixXd matIY0(const Eigen::VectorXd& Y) {
  int n = Y.size();
  Eigen::MatrixXd matIY(n,n);

  for (int k = 0; k < n; k++) {
    matIY.col(k) = (Y.array() >= Y[k]).cast<double>();
  }
  return matIY;
}
// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::export]]
// // each column of matIY[,i] is weight centered at X_i
// Eigen::SparseMatrix<double> matIYs(const Eigen::VectorXd& Y) {
//   int n = Y.size();
//   typedef Eigen::SparseMatrix<double> SparseMatrixType;
//   
//   SparseMatrixType matIY(n, n);
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       matIY.insert(i, j) = (Y(i) >= Y(j)) ? 1.0 : 0.0;
//     }
//   }
//   return matIY;
// }
//Eigen::MatrixXd matbw = matBw0(X,h)

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
// Eigen implementation of Sc function
Eigen::VectorXd Scp(double t,
          const Eigen::MatrixXd& wsc,
          const Eigen::VectorXd& Y, 
          const Eigen::VectorXd& delta) {
  Eigen::VectorXd IYw = (Y.array() <= t).cast<double>();
  Eigen::VectorXd deltaw = 1.0-delta.array();
  Eigen::VectorXd weight = IYw.array()*deltaw.array();
  Eigen::MatrixXd neww = wsc.array().colwise() * weight.array();
  Eigen::MatrixXd result = 1.0-neww.array();
  
  return result.colwise().prod();
}


// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::export]]
// // This is IPCW method
// Eigen::VectorXd IPCW_v(double h,
//                        const Eigen::VectorXd& Y, 
//                        const Eigen::VectorXd& X,
//                        const Eigen::VectorXd& delta,
//                        const Eigen::VectorXi& nn_index,
//                        const Eigen::VectorXi& index_death) {
//   int n = Y.size();
//   //calculate matBw: each columns is the weight centered at X_j
//   Eigen::MatrixXd matBw = matBw0(X, h);
//   //calculate matrix weight for Y>Y_i
//   Eigen::MatrixXd matIY = matIY0(Y);
//   //calculate the final weight matrix for Sc(t|X1),...,Sc(t|Xn)
//   Eigen::MatrixXd wsc = matBw.array()/(matIY.transpose()*matBw).array();
//   
//   Eigen::MatrixXd matsc(n,n);
//   for (int k = 0; k < n; ++k) {
//     matsc.row(k) = Scx(k,wsc,Y,delta);
//   }
//   Eigen::MatrixXd imatsc = 1.0/matsc.array();
//   
//   //calculate Q(Y_i) for Y_i is death.
//   Eigen::VectorXd Q_val = Eigen::VectorXd::Zero(n);
//   
//   for (int k = 0; k < index_death.size(); ++k) {
//     int index = index_death(k);
//     Q_val(index) = Qv(index,imatsc,Y,delta,nn_index);
//   }
//   double EQ = (Q_val.array()/matsc.diagonal().array()).mean();
//   double hat_IPCW = 6.0*EQ - 2.0;
//   // calculate BJ method
//   // we first calculate F(Y_i|X_i)
//   Eigen::MatrixXd B_Sc = matBw.array()/matsc.array();
//   Eigen::MatrixXd B_Sc_IY = B_Sc.array()*matIY.transpose().array();
//   Eigen::VectorXd Ft = delta.transpose()*B_Sc_IY; // this is all F(Y_i|X_i)
//   
//   // we then calculate Q(t)*F(t|X_i)
//   Eigen::MatrixXd B_Sc_iIY = B_Sc.array()*(1-matIY.array()).transpose().array();
//   Eigen::VectorXd Q_Ft = Q_val.transpose()*B_Sc_iIY; // this is all Q(t)*F(t|X_i)
//   // Q(t)*F(t|X_i)/(1-F(Y_i|X_i))
//   Eigen::VectorXd ratio = Q_Ft.array()/(1-Ft.array()).array();
//   //for BJ calculation
//   double sum_death = Q_val.dot(delta);
//   Eigen::VectorXd idelta = 1-delta.array();
//   double sum_censor = idelta.dot(ratio);
//   double mean_BJ = (sum_death+sum_censor)/static_cast<double>(n);
//   double BJ = 6.0*mean_BJ - 2.0;
//   Eigen::VectorXd result(2);
//   result << hat_IPCW,BJ;
//   return result;
// }
