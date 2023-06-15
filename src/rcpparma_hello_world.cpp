// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp"

#include "mixture.hpp"

using namespace Rcpp;

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
void compose_clusters_fast(const unsigned int nF, const unsigned int nH,
                           arma::mat& mu, arma::cube& sigma, arma::cube& sigma_inv,
                           const arma::cube& psi, const arma::mat& alpha) {
    // mu is nF + nH columns; column-major access faster
    // sigma is d by d by (nF + nH); slice-major access?
    // psi is d by d by nH; slice-major access?
    // alpha is nF by nH; column-major access faster
    for (unsigned int h = 0; h < nH; h++) {
        mu.col(nF + h) = mu.head_cols(nF) * alpha.col(h);
        
        sigma.slice(nF + h).zeros();
        for (unsigned int f = 0; f < nF; f++) sigma.slice(nF + h) += alpha(f, h) * alpha(f, h) * sigma.slice(f);
        sigma.slice(nF + h) += psi.slice(h);
        sigma.slice(nF + h) = (sigma.slice(nF + h) + sigma.slice(nF + h).t())/2;
        sigma_inv.slice(nF + h) = arma::inv_sympd(sigma.slice(nF + h));
    }
}

// [[Rcpp::export]]
arma::uvec draw_multiple(arma::mat probs) {
    unsigned int n = probs.n_rows;
    unsigned int g = probs.n_cols;
    
    arma::uvec proposed_assignment(n);
    for (unsigned int i = 0; i < n; i++) {
        double u = R::runif(0, 1);
        for (unsigned int j = 0; j < g; j++) {
            u -= probs(i,j);
            if (u < 0) {
                proposed_assignment(i) = j;
                break;
            }
        }
    }
    return proposed_assignment;
}

// [[Rcpp::export]]
void maximization_step_1_fast(const unsigned int nF, const unsigned int nH,
                              arma::mat& mu, arma::cube& sigma, arma::cube& sigma_inv,
                              arma::cube& psi, const arma::mat& alpha,
                              arma::vec& pi, const arma::mat& zhat,
                              const arma::mat& data, const String psi_type,
                              const String cov_type) {
    const unsigned int d = data.n_cols;
    unsigned int f, h, g;
    
    pi = arma::sum(zhat, 0).t();
    const arma::vec zhatf = pi;
    pi = arma::normalise(pi, 1);
    
    arma::field<arma::mat> ybar(nF, nH);
    for (h = 0; h < nH; h++) {
        arma::mat temp = (data.each_row() - mu.col(nF + h).t());
        for (f = 0; f < nF; f++) {
            ybar(f, h) = temp * (alpha(f, h) * sigma_inv.slice(nF + h) * sigma.slice(f));
            ybar(f, h).each_row() += mu.col(f).t();
        }
    }
    
    arma::field<arma::mat> Shfg(nH, nF, nF);
    for (h = 0; h < nH; h++) {
        for (f = 0; f < nF; f++) {
            for (g = 0; g < nF; g++) {
                Shfg(h, f, g) = (-(alpha(f, h) * alpha(g, h))) * sigma.slice(f) * sigma_inv.slice(nF + h) * sigma.slice(g);
            }
            Shfg(h, f, f) += sigma.slice(f);
        }
    }
    
    // mu update
    const double zhatftailnH = arma::accu(zhatf.tail(nH));
    for (f = 0; f < nF; f++) {
        mu.col(f) = data.t() * zhat.col(f);
        for (h = 0; h < nH; h++) mu.col(f) += ybar(f, h).t() * zhat.col(nF + h);
        mu.col(f) /= (zhatf(f) + zhatftailnH);
    }
    
    // sigma update
    const arma::mat sqrtzhat = arma::sqrt(zhat);
    // Function model_type("model.type");
    
    for (f = 0; f < nF; f++) {
        arma::mat temp = sqrtzhat.col(f) % (data.each_row() - mu.col(f).t()).each_col();
        sigma.slice(f) = temp.t() * temp;
        
        for (h = 0; h < nH; h++) {
            temp = sqrtzhat.col(nF + h) % (ybar(f,h).each_row() - mu.col(f).t()).each_col();
            sigma.slice(f) += zhatf(nF + h) * Shfg(h,f,f) + temp.t() * temp;
        }
        sigma.slice(f) /= (zhatf[f] + zhatftailnH);
        sigma.slice(f) = (sigma.slice(f) + sigma.slice(f).t())/2;
    }
    
    arma::cube factor_sigma = sigma.head_slices(nF);
    arma::vec factor_pi = pi.head(nF);
    
    arma::field<arma::cube> parsimonious = model_type_cpp(cov_type, factor_sigma, factor_pi);
    sigma.head_slices(nF) = parsimonious(0);
    sigma_inv.head_slices(nF) = parsimonious(1);
    
    const double min_psi_value = 1e-12;
    
    if (psi_type == "EV") {
        for (h = 0; h < nH; h++) {
            arma::mat tmp = data;
            for (f = 0; f < nF; f++) tmp -= alpha(f, h) * ybar(f, h);
            
            arma::mat temp = sqrtzhat.col(nF + h) % tmp.each_col();
            arma::vec psi_new = arma::sum(temp % temp, 0).t() / zhatf(nF + h);
            
            for (f = 0; f < nF; f++) for (g = 0; g < nF; g++) 
                psi_new += (alpha(f, h) * alpha(g, h)) * Shfg(h,f,g).diag();
            
            psi.slice(h).diag() = arma::clamp(psi_new, min_psi_value, INFINITY);
        }
    } else if (psi_type == "EE") { // diagonal, equal
        arma::mat psi_new(d, d, arma::fill::zeros);
        for (h = 0; h < nH; h++) {
            arma::mat tmp = data;
            for (f = 0; f < nF; f++) tmp -= alpha(f, h) * ybar(f, h);
            
            arma::mat temp = sqrtzhat.col(nF + h) % tmp.each_col();
            psi_new += (temp.t() * temp);
            
            for (f = 0; f < nF; f++) for (g = 0; g < nF; g++)
                psi_new += zhatf(nF + h) * (alpha(f, h) * alpha(g, h)) * Shfg(h, f, g);
        }
        psi_new /= zhatftailnH;
        psi_new = arma::diagmat(psi_new);
        psi.each_slice() = psi_new;
    } else if (psi_type == "IV") {
        for (h = 0; h < nH; h++) {
            arma::mat tmp = data;
            for (f = 0; f < nF; f++) tmp -= alpha(f, h) * ybar(f, h);
            
            arma::mat temp = sqrtzhat.col(nF + h) % tmp.each_col();
            arma::mat psi_new = (temp.t() * temp) / zhatf(nF + h);
            
            for (f = 0; f < nF; f++) for (g = 0; g < nF; g++) 
                psi_new += alpha(f, h) * alpha(g, h) * Shfg(h,f,g);
            
            psi.slice(h) = arma::eye(arma::size(psi_new)) * arma::trace(psi_new) / data.n_cols;
        }
    } else if (psi_type == "IE") {
        arma::mat psi_new(d, d, arma::fill::zeros);
        
        arma::mat tmp(size(data));
        for (h = 0; h < nH; h++) {
            tmp = data;
            for (f = 0; f < nF; f++) tmp -= alpha(f, h) * ybar(f, h);
            tmp.each_col() %= sqrtzhat.col(nF + h);
            
            psi_new += (tmp.t() * tmp);
            for (f = 0; f < nF; f++)
                for (g = 0; g < nF; g++)
                    psi_new += zhatf(nF + h) * alpha(f, h) * alpha(g, h) * Shfg(h,f,g);
        }
        
        psi_new /= zhatftailnH;
        psi_new = arma::eye(d, d) * arma::trace(psi_new) / d;
        psi.each_slice() = psi_new;
    } else if (psi_type == "C") {
        // Assume psi's are initialised correctly, don't update them.
    }
    return;
}

// [[Rcpp::export]]
void maximization_step_2_fast(const unsigned int nF, const unsigned int nH,
                              const arma::mat& mu, const arma::cube& sigma, const arma::cube& sigma_inv,
                              const arma::cube& psi, arma::mat& alpha,
                              const arma::vec& pi, const arma::mat& zhat,
                              const arma::mat& data, const arma::mat& Amat,
                              const arma::mat& bvec) {
    const unsigned int d = data.n_cols;
    const unsigned int n = data.n_rows;
    const arma::vec zhatf = arma::sum(zhat, 0).t();
    
    arma::field<arma::mat> ybar(nF, nH);
    for (unsigned int h = 0; h < nH; h++) {
        arma::mat temp = (data.each_row() - mu.col(nF + h).t());
        for (unsigned int f = 0; f < nF; f++) {
            ybar(f, h) = temp * (alpha(f, h) * sigma_inv.slice(nF + h) * sigma.slice(f));
            ybar(f, h).each_row() += mu.col(f).t();
        }
    }
    
    Function solve_QP("solve.QP");
    
    const arma::mat sqrtzhat = arma::sqrt(zhat);
    
    for (unsigned int h = 0; h < nH; h++) {
        // if (zhatf[nF + h] < 0.01) continue;
        
        arma::mat psiinv = arma::diagmat(1/psi.slice(h).diag());
        
        arma::mat Bc(nF, nF, arma::fill::zeros);
        for (unsigned int f = 0; f < nF; f++) for (unsigned int g = f; g < nF; g++) {
            Bc(f, g) = -alpha(f, h) * alpha(g, h) *
                arma::trace(psiinv * sigma.slice(f) * sigma_inv.slice(nF + h) * sigma.slice(g));
            if (f == g) {Bc(f, g) += arma::trace(psiinv * sigma.slice(f));}
            else Bc(g, f) = Bc(f, g);
        }
        arma::cube Pf(n, d, nF);
        for (unsigned int f = 0; f < nF; f++) Pf.slice(f) = ybar(f, h);
        Pf.each_slice() %= sqrtzhat.col(nF + h) * arma::sqrt(psiinv.diag()).t();
        Pf = arma::reshape(Pf, n * d, nF, 1);
        arma::mat Pf2 = Pf.slice(0).t() * Pf.slice(0) + zhatf(nF + h) * Bc;
        
        arma::mat qh_temp = data % (zhat.col(nF + h) * psiinv.diag().t());
        arma::vec qh(nF);
        for (unsigned int f = 0; f < nF; f++) qh(f) = arma::accu(ybar(f, h) % qh_temp);
        
        double max_val = std::max((arma::abs(Pf2)).max(), (arma::abs(qh)).max());
        
        Pf2 = Pf2 / max_val;
        qh = qh / max_val;
        
        Rcpp::List QP = solve_QP(Pf2, qh, Amat, bvec, 1);
        alpha.col(h) = as<arma::vec>(QP["solution"]);
    }
}

// [[Rcpp::export]]
arma::vec log_phi(arma::mat x, arma::vec mu, arma::mat sigma) {
    arma::mat x_cen = x.t();
    x_cen.each_col() -= mu;
    arma::solve(x_cen, arma::trimatl(arma::chol(sigma).t()), x_cen);
    x_cen = x_cen % x_cen;
    arma::vec const distval = arma::sum(x_cen, 0).t();    
    
    double const logdet = sum(arma::log(arma::eig_sym(sigma)));
    arma::vec const logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    return logretval;
}

// [[Rcpp::export]]
arma::vec log_phi_inv(arma::mat x, arma::vec mu, arma::mat sigma_inv) {
    arma::mat x_cen = (x.each_row() - mu.t()) * (arma::chol(sigma_inv).t());
    return -(x.n_cols * log2pi - arma::sum(arma::log(arma::eig_sym(sigma_inv))) + arma::sum(x_cen % x_cen, 1))/2;
}

// [[Rcpp::export]]
void expectation_step_fast(const unsigned int nF, const unsigned int nH,
                           const arma::mat& mu, const arma::cube& sigma, const arma::cube& sigma_inv,
                           const arma::vec& pi, arma::mat& zhat,
                           const arma::mat& data, const bool& hold_z,
                           const bool& stochastic, double& loglik) {
    
    arma::mat old_zhat = zhat;
    for (unsigned int g = 0; g < nF + nH; g++) {
        arma::mat x_cen = (data.each_row() - mu.col(g).t()) * (arma::chol(sigma_inv.slice(g)).t());
        zhat.col(g) = -(data.n_cols * log2pi - arma::sum(arma::log(arma::eig_sym(sigma_inv.slice(g)))) + arma::sum(x_cen % x_cen, 1))/2;
    }
    zhat.each_row() += log(pi.t());
    
    zhat = arma::exp(zhat);
    loglik = arma::accu(log(arma::sum(zhat, 1)));
    if (hold_z) {zhat = old_zhat;} else {zhat = arma::normalise(zhat, 1, 1);}
    
    return;
    
    // Ignore Stochastic Flag
    
    // if (stochastic) {
    //     arma::uvec assignments = draw_multiple(zhat);
    //     zhat.zeros();
    //     for (unsigned int n = 0; n < data.n_rows; n++) zhat(n, assignments(n)) = 1;
    // 
    //     // Check if each cluster is responsible for at least d observations
    //     arma::vec tabulated = arma::sum(zhat, 0).t();
    //     if (any(tabulated <= data.n_cols)) {
    //         zhat = old_zhat; // hold zhat as last iteration
    //     } // else keep zhat's the same
    // }
}

// [[Rcpp::export]]
NumericVector cpp_iteration_fast_progress(Rcpp::List state, 
                                          const unsigned int times,
                                          const unsigned int verbosity,
                                          const bool show_progress,
                                          const double AECM_epsilon = 1e-2) {
    
    // Unpack state list
    
    unsigned int nH = state["nH"];
    unsigned int nF = state["nF"];
    String psi_type = state["psi.type"];
    String cov_type = state["cov.type"];
    
    arma::mat mu = state["mu"];
    arma::cube sigma = state["sigma"];
    arma::cube sigma_inv(arma::size(sigma));
    for (unsigned int g = 0; g < nF + nH; g++) sigma_inv.slice(g) = arma::inv_sympd(sigma.slice(g));
    arma::mat alpha = state["alpha"];
    arma::cube psi = state["psi"];
    arma::mat zhat = state["zhat"];
    
    arma::vec pi = state["pi"];
    
    arma::mat data = state["data"];
    bool hold_z = state["hold.z"];
    bool stochastic = state["stochastic"];
    // bool stochastic = false;
    
    // Generate Amat and bvec here
    // These represent the simplicial constraint for the quadprog solver
    // Re-used often, so just make it once and re-use
    Rcpp::List convenience = state["convenience"];
    arma::mat Amat = convenience["Amat"];
    arma::vec bvec = convenience["bvec"];
    
    // Run iterations
    
    NumericVector ll(times);
    
    Progress p(times, show_progress);
    
    double ll_inf_old = -INFINITY;
    unsigned int i = 0;
    bool converged = false;
    
    if (verbosity == 0) {
        for (i = 0; i < times; i++) {
            maximization_step_1_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, psi_type, cov_type);
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            maximization_step_2_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, Amat, bvec);
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            p.increment();
            
            // Aitken's acceleration convergence criteria
            if (i > 2) {
                double at = (ll(i) - ll(i - 1)) / (ll(i - 1) - ll(i - 2));
                double ll_inf = ll(i - 1) + (ll(i) - ll(i - 1)) / (1 - at);
                if (std::abs(ll_inf - ll_inf_old) < AECM_epsilon) {
                    // Converged.
                    converged = true;
                    break;
                }
                ll_inf_old = ll_inf;
            }
            
        }
    } else if (verbosity == 1) {
        for (unsigned int i = 0; i < times; i++) {
            Rprintf("1");
            maximization_step_1_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, psi_type, cov_type);
            Rprintf("2");
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            Rprintf("3");
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            Rprintf("4");
            maximization_step_2_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, Amat, bvec);
            Rprintf("5");
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            Rprintf("6");
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            Rprintf("\b\b\b\b\b\b------\b\b\b\b\b\b");
            
            p.increment();
            // Aitken's acceleration convergence criteria
            if (i > 2) {
                double at = (ll(i) - ll(i - 1)) / (ll(i - 1) - ll(i - 2));
                double ll_inf = ll(i - 1) + (ll(i) - ll(i - 1)) / (1 - at);
                if (std::abs(ll_inf - ll_inf_old) < AECM_epsilon) {
                    // Converged.
                    converged = true;
                    break;
                }
                ll_inf_old = ll_inf;
            }
        }
    } else if (verbosity == 2) {
        arma::uvec subtiming(4, arma::fill::zeros);
        for (unsigned int i = 0; i < times; i++) {
            auto stage1 = std::chrono::high_resolution_clock::now();
            maximization_step_1_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, psi_type, cov_type);
            auto stage2 = std::chrono::high_resolution_clock::now();
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            auto stage3 = std::chrono::high_resolution_clock::now();
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            auto stage4 = std::chrono::high_resolution_clock::now();
            maximization_step_2_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha, pi, zhat, data, Amat, bvec);
            auto stage5 = std::chrono::high_resolution_clock::now();
            compose_clusters_fast(nF, nH, mu, sigma, sigma_inv, psi, alpha);
            auto stage6 = std::chrono::high_resolution_clock::now();
            expectation_step_fast(nF, nH, mu, sigma, sigma_inv, pi, zhat, data, hold_z, stochastic, ll(i));
            auto stage7 = std::chrono::high_resolution_clock::now();
            
            subtiming(0) += std::chrono::duration_cast<std::chrono::microseconds>(stage2 - stage1).count();
            subtiming(1) += std::chrono::duration_cast<std::chrono::microseconds>(stage3 - stage2).count();
            subtiming(2) += std::chrono::duration_cast<std::chrono::microseconds>(stage4 - stage3).count();
            subtiming(3) += std::chrono::duration_cast<std::chrono::microseconds>(stage5 - stage4).count();
            subtiming(1) += std::chrono::duration_cast<std::chrono::microseconds>(stage6 - stage5).count();
            subtiming(2) += std::chrono::duration_cast<std::chrono::microseconds>(stage7 - stage6).count();
            
            p.increment();
            // Aitken's acceleration convergence criteria
            if (i > 2) {
                double at = (ll(i) - ll(i - 1)) / (ll(i - 1) - ll(i - 2));
                double ll_inf = ll(i - 1) + (ll(i) - ll(i - 1)) / (1 - at);
                if (std::abs(ll_inf - ll_inf_old) < AECM_epsilon) {
                    // Converged.
                    converged = true;
                    break;
                }
                ll_inf_old = ll_inf;
            }
        }
        
        subtiming.print("Subtimings");
    }
    
    // Repack state list
    
    state["mu"] = mu;
    state["sigma"] = sigma;
    state["alpha"] = alpha;
    state["psi"] = psi;
    state["zhat"] = zhat;
    state["pi"] = pi;
    state["loglik"] = ll(i - 1);
    state["converged"] = converged;
    
    ll = head(ll, i);
    
    return ll;
}
