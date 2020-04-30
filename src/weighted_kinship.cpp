#include <Rcpp.h>
#include <RcppParallel.h> 
#include <iostream>
#include "gaston/matrix4.h"
using namespace Rcpp;
using namespace RcppParallel;


// ************* [on ne symmétrise pas] ***********

struct paraKin_weighted : public Worker {
  // input and others
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol;
  const double * mu;
  const double * w;
  const size_t sizeK;

  // output
  double * K;
  
  // constructeurs
  paraKin_weighted(uint8_t ** data, const size_t ncol, const size_t true_ncol, const double * mu, const double * w) : 
        data(data), ncol(ncol), true_ncol(true_ncol), mu(mu), w(w), sizeK((4*true_ncol)*(4*true_ncol+1)/2) { 
          K = new double[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
        }
  paraKin_weighted(paraKin_weighted & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), mu(Q.mu), w(Q.w), sizeK(Q.sizeK) {
          K = new double[sizeK];  // K is padded to a multiple of 4...
          std::fill(K, K+sizeK, 0);
        }

  // destructeur
  ~paraKin_weighted() { 
          delete [] K; 
  }

  // worker !
  void operator()(size_t beg, size_t end) {
    double H[4];
    double h0[32];
    double h1[32];
    H[3] = 0;

    for(size_t i = beg; i < end; i++) {
      double w_ = (double) w[i]; 
      if(w_ == 0) continue;
      double mu_ = (double) mu[i];
      double v0 = -mu_;
      double v1 = (1-mu_);
      double v2 = (2-mu_);
      H[0] = v0;
      H[1] = v1;
      H[2] = v2;
      
      size_t k = 0;
      for(size_t a1 = 0; a1 < 4; a1++) {
        for(size_t a2 = 0; a2 < 4; a2++) {
          h0[k++] = H[a2]*w_;
          h0[k++] = H[a1]*w_;
        }
      }

      const uint8_t * dd = data[i];
      
      k = 0;
      for(size_t j1 = 0; j1 < true_ncol; j1++) {
        uint8_t x1 = dd[j1];
        for(unsigned int ss1 = 0; (ss1 < 4); ss1++) {
          for(size_t a = 0; a < 32; a++) h1[a] = H[x1&3]*h0[a];
          for(size_t j2 = 0; j2 < j1; j2++) {
            uint8_t x2 = dd[j2];
            K[k++] += h1[ (x2&15)<<1 ];
            K[k++] += h1[ ((x2&15)<<1)+1 ];
            K[k++] += h1[ (x2>>4)<<1 ];
            K[k++] += h1[ ((x2>>4)<<1)+1 ];
          }
          size_t j2 = j1;
          uint8_t x2 = dd[j2];
          for(unsigned int ss2 = 0; ss2 <= ss1; ss2++) {
            K[k++] += H[x1&3]*H[x2&3]*w_; 
            x2 >>= 2;
          } 
          x1 >>= 2; 
        }
      } 
    }
  }

  // recoller
  void join(const paraKin_weighted & Q) {
    std::transform(K, K + sizeK, Q.K, K, std::plus<double>());
    // autrement dit : K += Q.K;
  }

};


// [[Rcpp::export]]
NumericMatrix weighted_Kinship_w(XPtr<matrix4> p_A, const std::vector<double> & mu, const std::vector<double> & w, LogicalVector snps, int chunk) {
  int nb_snps = sum(snps);

  if(snps.length() != p_A->nrow || mu.size() != nb_snps || w.size() != nb_snps) 
    stop("Dimensions mismatch");

  uint8_t ** data = new uint8_t * [nb_snps];
  size_t k = 0;
  for(size_t i = 0; i < p_A->nrow; i++) {
    if(snps[i]) data[k++] = p_A->data[i];
  }

  paraKin_weighted X(data, p_A->ncol, p_A->true_ncol, &mu[0], &w[0]);
  parallelReduce(0, nb_snps, X, chunk);

  delete [] data;

  NumericMatrix Y(p_A->ncol,p_A->ncol);
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(j,i) = (double) X.K[k++];
    }
  }

  // symetriser
  k = 0;
  for(size_t i = 0; i < p_A->ncol; i++) {
    for(size_t j = 0; j <= i; j++) {
      Y(i,j) = X.K[k++]; // ou Y(j,i);
    }
  }
  
  return Y;
}


