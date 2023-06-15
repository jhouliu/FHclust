initialization = function(data, num.F, num.H, cov.type = "VVV", psi.type = "EV",
                          init.assign = NULL, psi.C = 1e-10, subsample = 0.5) {
  # Input Parameter Sanity Check
  if (!Rfast::is_integer(num.F)) stop("num.F is not an integer!")
  if (!Rfast::is_integer(num.H)) stop("num.H is not an integer!")
  if (num.F < 2) stop("num.F must be at least 2!")
  if (num.H < 1) stop("num.H must be at least 1!")
  psi.type = toupper(psi.type)
  if (length(psi.type) != 1) stop("Invalid psi.type!")
  if (!(psi.type %in% c("EV", "EE", "IV", "IE", "C"))) stop("Unrecognized psi.type!")
  
  # Initialize state of fit
  state = list(data = data, nF = num.F, nH = num.H,
               zhat = matrix(1/(num.F + num.H), nrow = nrow(data), ncol = num.F + num.H),
               pi = rep(1/(num.F + num.H), num.F + num.H),
               psi.type = psi.type,
               cov.type = cov.type,
               stochastic = FALSE,
               hold.z = FALSE,
               converged = FALSE,
               convenience = list(Amat = cbind(1L, diag(num.F)),
                                  bvec = c(1L, rep(0, num.F)),
                                  magic.mat = suppressWarnings(matrix(c(rep(1, ncol(data)), rep(0, ncol(data) * num.F)), ncol(data) * num.F, num.F))
               ) 
  )
  
  half.data = data[sample.int(nrow(data), size = nrow(data) * subsample),]
  
  trial.G = num.F + num.H
  
  repeat {
    mcfit = Mclust(half.data, trial.G, cov.type, verbose = FALSE)
    if (is.null(mcfit)) {
      # Release cov.type restriction
      mcfit = Mclust(half.data, trial.G, verbose = FALSE)
      mcll  = mclustLoglik(mc$BIC)[1,]
      # mcbic = mc$BIC[1,]
      # mcdf  = sapply(names(mcbic), nMclustParams, d = ncol(data), G = trial.G)
      # mcind = which(names(mcbic) == cov.type)
      # mcbic = mcbic[1:mcind]
      available.cov.type = names(which.max(mcll))
      if (is.null(available.cov.type)) stop("Failed to initialize.")
      mcfit = Mclust(half.data, trial.G, available.cov.type, verbose = FALSE)
    }
      
    ext.pts = if (ncol(data) + 1 > trial.G) 1:trial.G else 
      ext.pts = (unique(as.vector(convhulln(t(mcfit$parameters$mean)))))
    if (length(ext.pts) >= num.F) break
    trial.G = trial.G + 1L
    if (trial.G > 10 * (num.F + num.H)) stop("Failed to initialize; couldn't get num.F extremal points.")
  }
  
  if (length(ext.pts) > num.F) {
    trial.indices = combn(length(ext.pts), num.F)
    discrepancies = apply(trial.indices, 2, function(indices) {
      M = mcfit$parameters$mean[,ext.pts[indices]]
      Dmat = forcePD(crossprod(M))
      maxval = max(abs(Dmat))
      Dmat = Dmat / maxval
      devs = sapply(ext.pts[-indices], function(x) {
        dvec = colsums(mcfit$parameters$mean[,x] * M)
        dvec = dvec / maxval
        dev = solve.QP(Dmat, dvec, 
                       cbind(1, diag(num.F)), c(1, rep(0, num.F)), 1)$value * 2 * maxval +
          sum(mcfit$parameters$mean[,ext.pts[x]]^2)
        return(dev)
      })
      return(sum(devs))
    })
    num.F.select = ext.pts[trial.indices[,which.min(discrepancies)]]
  } else {
    # Just enough exterior points
    num.F.select = ext.pts
  }
  
  # Get num.H from the remaining clusters
  num.H.select = setdiff(1:trial.G, num.F.select)
  # Grab clusters with most responsibility from non-prototype clusters
  num.H.select = num.H.select[head(order(mcfit$parameters$pro[num.H.select]), num.H)]
  
  # Estimate alphas for these compositional clusters
  alpha.list = lapply(num.H.select, function(c) {
    M = mcfit$parameters$mean[, num.F.select]
    Dmat = forcePD(crossprod(M))
    dvec = colsums(mcfit$parameters$mean[,c] * M)
    maxval = max(abs(Dmat), abs(dvec))
    Dmat = Dmat / maxval
    dvec = dvec / maxval
    alphas = solve.QP(Dmat, dvec, 
                      cbind(1, diag(num.F)), c(1, rep(0, num.F)), 1)$solution
    return(alphas)
  })
  
  # Put together state object
  state$mu                                  = matrix(0, nrow = ncol(data), ncol = num.F + num.H)
  state$mu[, 1:num.F]                       = mcfit$parameters$mean[, num.F.select]
  state$sigma                               = array(0, dim = c(ncol(data), ncol(data), num.F + num.H))
  state$sigma[, , 1:num.F]                  = mcfit$parameters$variance$sigma[, , num.F.select] * 1
  state$sigma_inv                           = state$sigma
  for (p in 1:num.F) state$sigma_inv[, , p] = solve(state$sigma_inv[, , p]); rm(p)
  state$psi                                 = array(0, dim = c(ncol(data), ncol(data), num.H))
  
  cov.diag = diag(cov(half.data))
  
  state$psi = switch(psi.type,
                     EV = array(diag(cov.diag), dim = c(ncol(data), ncol(data), num.H)),
                     EE = array(diag(cov.diag), dim = c(ncol(data), ncol(data), num.H)),
                     IV = array(diag(mean(cov.diag), ncol(data)), dim = c(ncol(data), ncol(data), num.H)),
                     IE = array(diag(mean(cov.diag), ncol(data)), dim = c(ncol(data), ncol(data), num.H)),
                     C = array(psi.C * diag(ncol(data)), dim = c(ncol(data), ncol(data), num.H)))
  # state$psi             = array(psi.II * diag(ncol(data)), dim = c(ncol(data), ncol(data), num.H))
  state$alpha           = do.call(cbind, alpha.list)
  
  compose_clusters_fast(num.F, num.H, state$mu, state$sigma, state$sigma_inv, state$psi, state$alpha)
  
  loglik = 0. # Dummy placeholder
  expectation_step_fast(num.F, num.H, state$mu, state$sigma, state$sigma_inv, state$pi, state$zhat,
                        state$data, state$hold.z, state$stochastic, loglik)
  state$loglik = 0
  
  # Apply Parsimonious Covariance Matrices
  parsimonious_cov = model_type_cpp(cov.type, state$sigma[, , 1:num.F], state$pi[1:num.F])
  state$sigma[, , 1:num.F] = parsimonious_cov[[1]]
  state$sigma_inv[,,1:num.F] = parsimonious_cov[[2]]
  compose_clusters_fast(num.F, num.H, state$mu, state$sigma, state$sigma_inv, state$psi, state$alpha)
  expectation_step_fast(num.F, num.H, state$mu, state$sigma, state$sigma_inv, state$pi, state$zhat,
                        state$data, state$hold.z, state$stochastic, loglik)
  
  return(state)
}