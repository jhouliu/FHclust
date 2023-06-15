// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::field<arma::cube> getEkOk_cpp(const arma::cube Sk, const arma::vec ng, const unsigned int G) {
  arma::cube Ok(arma::size(Sk));
  arma::cube EWk = Sk;
  
  for (unsigned int g = 0; g < G; g++) {
    arma::mat Wk = Sk.slice(g) * ng(g);
    
    arma::vec eigvals; // never used in this function
    arma::mat eigvecs;
    arma::eig_sym(eigvals, eigvecs, Wk);
    
    EWk.slice(g) = arma::fliplr(eigvecs); // order of eigenvalues is reversed
    Ok.slice(g) = EWk.slice(g).t() * Wk * EWk.slice(g);
  }
  
  arma::field<arma::cube> result(2);
  result(0) = Ok;
  result(1) = EWk;
  return result;
}

// [[Rcpp::export]]
arma::mat getA_cpp(const arma::cube Ok, const arma::vec lam, const unsigned int G, const unsigned int d) {
  arma::mat A(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) A += Ok.slice(g) / lam(g);
  arma::vec Adiag = A.diag();
  A = arma::diagmat(Adiag / pow(arma::prod(Adiag), 1.0/d));
  return A;
}

// [[Rcpp::export]]
arma::mat sumSk_wt_cpp(const arma::cube Sk, const arma::vec wt, const unsigned int d, const unsigned int G) {
  arma::mat W(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) W += Sk.slice(g) * wt(g);
  return W;
}

// [[Rcpp::export]]
arma::mat newD3_MM_cpp(const arma::mat D, const unsigned int d, const unsigned int G, const arma::cube Wk, const arma::mat Ak, const unsigned int tmax = 100) {
  arma::mat z(d, d, arma::fill::zeros);
  // double lambda = 0;
  
  for (unsigned int g = 0; g < G; g++) {
    arma::vec eigvals;
    arma::mat eigvecs; // never used in this function
    arma::eig_sym(eigvals, eigvecs, Wk.slice(g));
    
    double lambdak = arma::max(eigvals);
    z += arma::diagmat(1.0 / Ak.col(g)) * D.t() * Wk.slice(g) -
      lambdak * (arma::diagmat(1.0 / Ak.col(g)) * D.t());
  }
  
  arma::mat svd_u, svd_v; arma::vec svd_d; // pieces of z1
  svd(svd_u, svd_d, svd_v, z);
  
  arma::mat Xk1 = svd_v * svd_u.t();
  return Xk1;
}

// [[Rcpp::export]]
arma::mat newD4_MM_cpp(const arma::mat D, const unsigned int d, const unsigned int G, const arma::cube Wk, const arma::mat Ak, const unsigned int tmax = 100) {
  arma::mat z(d, d, arma::fill::zeros);
  // double lambda = 0;
  
  for (unsigned int g = 0; g < G; g++) {
    double lambdak = arma::max(1 / Ak.col(g));
    z += Wk.slice(g) * D * arma::diagmat(1 / Ak.col(g)) -
      lambdak * (Wk.slice(g) * D);
  }
  
  arma::mat svd_u, svd_v; arma::vec svd_d; // pieces of z1
  svd(svd_u, svd_d, svd_v, z);
  arma::mat Xk1 = (svd_v * svd_u.t()).t(); // transpose rolled into this line
  return Xk1;
}

// [[Rcpp::export]]
arma::mat newD_cpp(const arma::mat D, const unsigned int d, const unsigned int G, const arma::cube Wk, const arma::mat Ak, const unsigned int tmax = 100) {
  arma::mat D6 = D;
  D6 = newD3_MM_cpp(D6, d, G, Wk, Ak, tmax); // Pass-through tmax
  D6 = newD4_MM_cpp(D6, d, G, Wk, Ak, tmax); // Pass-through tmax
  return D6;
}

// [[Rcpp::export]]
double testval_cpp(const arma::cube Wk, const arma::mat Ak, const arma::mat D, const unsigned int G) {
  arma::vec z(G);
  for (unsigned int g = 0; g < G; g++) z(g) = arma::trace(D.t() * Wk.slice(g) * D * arma::diagmat(1/Ak.col(g)));
  double result = arma::accu(z);
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEEE_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::mat W = sumSk_wt_cpp(Sk, ng, d, G) / arma::accu(ng);
  
  // Remove dead code
  // arma::cube val(d,d,G);
  // val.each_slice = W;
  
  double logdetW = log(arma::det(W));
  arma::mat invW = arma::inv(W);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  sigma.each_slice() = W;
  invSigma.each_slice() = invW;
  logdet.fill(logdetW);
  // Replaced by the slightly faster code above.
  // 
  // for (unsigned int g = 0; g < G; g++) {
  //   sigma.slice(g) = W;
  //   invSigma.slice(g) = invW;
  //   logdet.slice(g) = logdetW;
  // }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEEV_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-12, const unsigned int max_iter = 100) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube EWk = Sk;
  arma::mat A(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) {
    arma::mat Wk = Sk.slice(g) * ng(g);
    
    arma::vec eigvals; // never used in this function
    arma::mat eigvecs;
    arma::eig_sym(eigvals, eigvecs, Wk);
    
    EWk.slice(g) = arma::fliplr(eigvecs); // order of eigenvalues is reversed
    A += EWk.slice(g).t() * Wk * EWk.slice(g);
  }
  
  double lam = pow(arma::prod(A.diag()), 1.0 / d);
  A /= lam;
  lam /= arma::accu(ng);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam * (EWk.slice(g) * A * EWk.slice(g).t());
    invSigma.slice(g) = 1/lam * (EWk.slice(g) * arma::diagmat(1.0 / A.diag()) * EWk.slice(g).t());
    logdet.slice(g) = d * log(lam);
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}


// [[Rcpp::export]]
arma::field<arma::cube> msVEV_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-14, const unsigned int max_iter = 100) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::field<arma::cube> temp = getEkOk_cpp(Sk, ng, G);
  arma::cube Ok = temp(0);
  arma::cube EWk = temp(1);
  arma::vec lam(G); for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Ok.slice(g)) / (ng(g) * d);
  arma::mat A = getA_cpp(Ok, lam, G, d);
  arma::mat invA(size(A), arma::fill::zeros);
  invA.diag() = 1.0 / A.diag();
  for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Ok.slice(g) % invA) / (ng(g) * d);
  
  double conv1 = d * arma::accu(ng % (1 + log(lam)));
  double conv2 = INFINITY;
  
  unsigned int count = 1;
  while ((conv2 - conv1) / conv1 > eplison && count < max_iter) {
    A = getA_cpp(Ok, lam, G, d);
    invA.diag() = 1.0 / A.diag();
    for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Ok.slice(g) % invA) / (ng(g) * d);
    conv2 = conv1;
    conv1 = d * arma::accu(ng % (1.0 + log(lam)));
    count += 1;
  }
  
  // dead code in original msVEV function
  // arma::cube val(d, d, G, arma::fill::zeros);
  // for (unsigned int g = 0; g < G; g++) val.slice(g) = lam(g) * (EWk.slice(g) * A * EWk.slice(g).t());
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam(g) * (EWk.slice(g) * A * EWk.slice(g).t());
    invSigma.slice(g) = (1.0 / lam(g)) * (EWk.slice(g) * invA * EWk.slice(g).t());
    logdet.slice(g) = d * log(lam(g));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVVV_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = Sk.slice(g);
    invSigma.slice(g) = arma::inv(Sk.slice(g));
    logdet.slice(g) = log(arma::det(Sk.slice(g)));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEEI_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::mat W = sumSk_wt_cpp(Sk, ng, d, G) / arma::accu(ng);
  arma::mat B = arma::diagmat(W);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  sigma.each_slice() = B;
  invSigma.each_slice() = arma::diagmat(1 / B.diag());
  logdet.fill(arma::accu(log(B.diag())));
  // Replaced by the slightly faster code above.
  // 
  // arma::mat invB = arma::diagmat(1 / W.diag());
  // for (unsigned int g = 0; g < G; g++) {
  //   sigma.slice(g) = B;
  //   invSigma.slice(g) = invB;
  //   logdet.slice(g) = arma::accu(log(B.diag()));
  // }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVEI_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-20, const unsigned int max_iter = 100) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::vec lam(G); for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Sk.slice(g)) / d;
  arma::mat W = sumSk_wt_cpp(Sk, ng / lam, d, G);
  arma::mat B = arma::diagmat(W.diag() / pow(arma::prod(W.diag()), 1.0/d));
  arma::mat invB = arma::diagmat(1.0 / B.diag());
  for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Sk.slice(g) % invB);
  
  double conv1 = d * arma::accu(ng % (1 + log(lam)));
  double conv2 = INFINITY;
  
  unsigned int count = 1;
  while (abs(conv2 - conv1) > eplison && count < max_iter) {
    W = sumSk_wt_cpp(Sk, ng / lam, d, G);
    B = arma::diagmat(W.diag() / pow(arma::prod(W.diag()), 1.0/d));
    invB = arma::diagmat(1.0 / B.diag());
    for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Sk.slice(g) % invB) / d;
    conv2 = conv1;
    conv1 = d * arma::accu(ng % (1.0 + log(lam)));
    count += 1;
  }
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam(g) * B;
    invSigma.slice(g) = 1.0 / lam(g) * invB;
    logdet.slice(g) = d * log(lam(g));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEVI_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::mat Bk(d, G);
  arma::vec lam(G);
  for (unsigned int g = 0; g < G; g++) {
    Bk.col(g) = Sk.slice(g).diag() * ng(g);
    lam(g) = pow(arma::prod(Bk.col(g)), 1.0 / d);
    Bk.col(g) /= lam(g);
  }
  // Bk.each_row() /= lam.t();
  
  double lam2 = arma::accu(lam) / arma::accu(ng);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam2 * arma::diagmat(Bk.col(g));
    invSigma.slice(g) = 1.0 / lam2 * arma::diagmat(1.0 / Bk.col(g));
    logdet.slice(g) = d * log(lam2);
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVVI_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = arma::diagmat(Sk.slice(g));
    invSigma.slice(g) = arma::diagmat(1 / Sk.slice(g));
    logdet.slice(g) = arma::accu(log(Sk.slice(g).diag()));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEII_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::mat W = sumSk_wt_cpp(Sk, ng, d, G) / arma::accu(ng);
  double lam = arma::trace(W) / (arma::accu(ng) * d);
  
  arma::cube sigma(d,d,G,arma::fill::zeros);
  arma::cube invSigma(d,d,G,arma::fill::zeros);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g).diag().fill(lam);
    invSigma.slice(g).diag().fill(1.0 / lam);
  }
  logdet.fill(d * log(lam));
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVII_cpp(const arma::cube Sk, const arma::vec ng) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube sigma(d,d,G,arma::fill::zeros);
  arma::cube invSigma(d,d,G,arma::fill::zeros);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    double sumdiagSkg = arma::accu(Sk.slice(g).diag());
    sigma.slice(g).diag().fill(sumdiagSkg / d);
    invSigma.slice(g).diag().fill(d / sumdiagSkg);
    logdet.slice(g) = d * log(sumdiagSkg) - d * log(d);
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVVE_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-20, const unsigned int max_iter = 10) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube Wk = Sk;
  arma::mat W(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) {
    Wk.slice(g) *= ng(g);
    W += Wk.slice(g);
  }
  
  // Assume hard-coded D0 = NULL
  arma::mat D = arma::eye(d, d);
  
  arma::mat Ak(d, G);
  for (unsigned int g = 0; g < G; g++) {
    arma::mat temp = D.t() * Wk.slice(g) * D; // Need to save this result first
    Ak.col(g) = temp.diag();
    Ak.col(g) /= pow(arma::prod(Ak.col(g)), 1.0 / d);
  }
  D = newD_cpp(D, d, G, Wk, Ak);
  
  double conv1 = testval_cpp(Wk, Ak, D, G);
  double conv2 = INFINITY;
  
  unsigned int count = 1;
  while ((conv2 - conv1)/abs(conv1) > eplison && count < max_iter) {
    for (unsigned int g = 0; g < G; g++) {
      arma::mat temp = D.t() * Wk.slice(g) * D; 
      Ak.col(g) = temp.diag();
      Ak.col(g) /= pow(arma::prod(Ak.col(g)), 1.0 / d);
    }
    D = newD_cpp(D, d, G, Wk, Ak);
    
    conv2 = conv1;
    conv1 = testval_cpp(Wk, Ak, D, G);
    count += 1;
  }
  
  arma::vec lam(G);
  for (unsigned int g = 0; g < G; g++) {
    lam(g) = arma::trace(D * arma::diagmat(1 / Ak.col(g)) * D.t() * Sk.slice(g)) / d;
  }
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = D * arma::diagmat(lam(g) * Ak.col(g)) * D.t();
    invSigma.slice(g) = D * arma::diagmat(1 / lam(g) * (1 / Ak.col(g))) * D.t();
    logdet.slice(g) = d * log(lam(g));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEVE_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-20, const unsigned int max_iter = 10) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube Wk = Sk;
  arma::mat W(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) {
    Wk.slice(g) *= ng(g);
    W += Wk.slice(g);
  }
  
  // Assume hard-coded D0 = NULL
  arma::mat D = arma::eye(d, d);
  
  arma::mat Ak(d, G);
  for (unsigned int g = 0; g < G; g++) {
    arma::mat temp = D.t() * Wk.slice(g) * D; // Need to save this result first
    Ak.col(g) = temp.diag();
    Ak.col(g) /= pow(arma::prod(Ak.col(g)), 1.0 / d);
  }
  D = newD_cpp(D, d, G, Wk, Ak);
  
  double conv1 = testval_cpp(Wk, Ak, D, G);
  double conv2 = INFINITY;
  
  unsigned int count = 1;
  while ((conv2 - conv1)/abs(conv1) > eplison && count < max_iter) {
    D = newD_cpp(D, d, G, Wk, Ak); // Unclear why this moved compared to VVE
    for (unsigned int g = 0; g < G; g++) {
      arma::mat temp = D.t() * Wk.slice(g) * D; 
      Ak.col(g) = temp.diag();
      Ak.col(g) /= pow(arma::prod(Ak.col(g)), 1.0 / d);
    }
    
    conv2 = conv1;
    conv1 = testval_cpp(Wk, Ak, D, G);
    count += 1;
  }
  
  double lam = 0;
  for (unsigned int g = 0; g < G; g++) {
    lam += arma::trace(D * arma::diagmat(1/Ak.col(g)) * D.t() * Wk.slice(g));
  }
  lam /= (arma::accu(ng) * d);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = D * arma::diagmat(lam * Ak.col(g)) * D.t();
    invSigma.slice(g) = D * arma::diagmat(1 / lam * (1 / Ak.col(g))) * D.t();
    logdet.slice(g) = d * log(lam);
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msVEE_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-14, const unsigned int max_iter = 100) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube Wk = Sk;
  arma::mat W(d, d, arma::fill::zeros);
  for (unsigned int g = 0; g < G; g++) {
    Wk.slice(g) *= ng(g);
    W += Wk.slice(g);
  }
  
  arma::mat C = W / pow(arma::det(W), 1.0 / d);
  arma::mat invC = arma::inv(C);
  arma::vec lam(G); for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Sk.slice(g) * invC) / d;
  
  double val1 = 0;
  for (unsigned int g = 0; g < G; g++) val1 += arma::trace(Wk.slice(g) % invC)/lam(g);
  val1 += d * arma::accu(ng % lam);
  
  double conv1 = val1;
  double conv2 = INFINITY;
  
  unsigned int count = 1;
  while ((conv2 - conv1)/abs(conv1) > eplison && count < max_iter) {
    C = sumSk_wt_cpp(Wk, 1.0 / lam, d, G);
    C = C / pow(arma::det(C), 1.0 / d);
    invC = arma::inv(C);
    
    for (unsigned int g = 0; g < G; g++) lam(g) = arma::trace(Sk.slice(g) * invC) / d;
    
    double val1 = 0;
    for (unsigned int g = 0; g < G; g++) val1 += arma::trace(Wk.slice(g) % invC)/lam(g);
    val1 += d * arma::accu(ng % lam);
    
    conv2 = conv1;
    conv1 = val1;
    count += 1;
  }
  
  invC = arma::inv(C);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam(g) * C;
    invSigma.slice(g) = 1.0/lam(g) * invC;
    logdet.slice(g) = d * log(lam(g));
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> msEVV_cpp(const arma::cube Sk, const arma::vec ng, const double eplison = 1e-12, const unsigned int max_iter = 100) {
  unsigned int d = Sk.n_rows;
  unsigned int G = Sk.n_slices;
  
  arma::cube Wk = Sk;
  arma::vec lam(G);
  
  for (unsigned int g = 0; g < G; g++) {
    Wk.slice(g) *= ng(g);
    lam(g) = pow(arma::det(Wk.slice(g)), 1.0/d);
  }
  
  arma::cube Ck = Wk;
  for (unsigned int g = 0; g < G; g++) Ck.slice(g) /= lam(g);
  double lam2 = arma::accu(lam);
  
  arma::cube sigma(d,d,G);
  arma::cube invSigma(d,d,G);
  arma::cube logdet(1,1,G);
  
  for (unsigned int g = 0; g < G; g++) {
    sigma.slice(g) = lam2 * Ck.slice(g);
    invSigma.slice(g) = 1.0/lam2 * arma::inv(Ck.slice(g));
    logdet.slice(g) = d * log(lam2);
  }
  
  arma::field<arma::cube> result(3);
  result(0) = sigma;
  result(1) = invSigma;
  result(2) = logdet;
  return result;
}

// [[Rcpp::export]]
arma::field<arma::cube> model_type_cpp (const Rcpp::String modelname, const arma::cube Sk, 
                                        arma::vec ng, const double mtol = 1e-18, const unsigned int mmax = 10) {
  arma::field<arma::cube> val;
  ng = arma::normalise(ng, 1);
  
  if (modelname == "EII") val = msEII_cpp(Sk, ng);
  else if (modelname == "VII") val = msVII_cpp(Sk, ng);
  else if (modelname == "EEI") val = msEEI_cpp(Sk, ng);
  else if (modelname == "VEI") val = msVEI_cpp(Sk, ng, mtol, mmax);
  else if (modelname == "EVI") val = msEVI_cpp(Sk, ng);
  else if (modelname == "VVI") val = msVVI_cpp(Sk, ng);
  else if (modelname == "EEE") val = msEEE_cpp(Sk, ng);
  else if (modelname == "EEV") val = msEEV_cpp(Sk, ng);
  else if (modelname == "VEV") val = msVEV_cpp(Sk, ng, mtol, mmax);
  else if (modelname == "VVV") val = msVVV_cpp(Sk, ng);
  else if (modelname == "EVE") val = msEVE_cpp(Sk, ng, mtol, mmax);
  else if (modelname == "VVE") val = msVVE_cpp(Sk, ng, mtol, mmax);
  else if (modelname == "VEE") val = msVEE_cpp(Sk, ng, mtol, mmax);
  else if (modelname == "EVV") val = msEVV_cpp(Sk, ng);
  else stop("modelname or covtype is not correctly defined");
  
  return val;
}
