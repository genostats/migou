// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
typedef Eigen::Map<Eigen::VectorXd> VEC;
typedef Eigen::Map<Eigen::MatrixXd> MAT;

// [[Rcpp::export]]
double min_var_eigen(NumericVector x_, NumericMatrix a_, int N = 100, double eps = 1e-6) {
  VEC x(as<VEC>(x_));
  MAT A(as<MAT>(a_));
  int n = A.rows();

  VectorXd g = A*x;
  double mean_g = g.mean();
  for(int i = 0; i < n; i++) g(i) -= mean_g;

  VectorXd y = g;
  VectorXd Ay = A*y;

  double yAy = y.dot(Ay);
  double alpha = x.dot(Ay) / yAy;
  x.noalias() -= alpha*y;

  for(int i = 0; i < N; i ++) {
    g.noalias() = A*x;
    double mean_g = g.mean();
    for(int i = 0; i < n; i++) g(i) -= mean_g;

    y = g - (g.dot(Ay)/yAy)*y;
    if( y.cwiseAbs().maxCoeff() < eps ) break;
    Ay.noalias() = A*y;
    yAy = y.dot(Ay);
    alpha = x.dot(Ay) / yAy;
    x.noalias() -= alpha*y;
  }
  return x.dot(A*x);
}
