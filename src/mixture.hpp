#ifndef PACKAGENAME_ADD_H
#define PACKAGENAME_ADD_H

// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

arma::field<arma::cube> model_type_cpp (const Rcpp::String modelname, const arma::cube Sk, const arma::vec ng, const double mtol = 1e-18, const unsigned int mmax = 10) ;

#endif
