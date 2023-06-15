summary_state = function(state, true.classes = NULL) {
  if (is.null(state)) {
    return(c(
      ll = NA,
      df = NA,
      bic = NA,
      ari = NA,
      min.eig = NA,
      min.z = NA,
      min.assign = NA,
      total.iters = NA,
      converged = NA
    ))
  }
  
  d = nrow(state$mu)
  df = state$nF * d + ncovpar(state$cov.type, d, state$nF) + 
    state$nH * (state$nF - 1) + 
    switch(state$psi.type,
           EV = state$nH * d,
           EE = d,
           IV = state$nH,
           IE = 1,
           C = 0) + 
    length(state$pi) - 1
  bic = log(nrow(state$zhat)) * df - 2 * state$loglik
  ari = if (is.null(true.classes)) NA else mclust::adjustedRandIndex(rowMaxs(state$zhat), true.classes)
  min.eig = min(apply(state$sigma, 3, function(x) eigen(x, only.values = TRUE)$values))
  return(c(
    ll = state$loglik,
    df = df,
    bic = bic,
    ari = ari,
    min.eig = min.eig,
    min.z = min(colSums(state$zhat)),
    min.assign = min(table(factor(max.col(state$zhat), 1:(state$nF + state$nH)))),
    total.iters = length(state$ll),
    converged = state$converged
  ))
}

summary.FHclust = function(state, true.classes = NULL, ...) {
  nF = state$nF
  nH = state$nH
  
  d = nrow(state$mu)
  df = state$nF * d + ncovpar(state$cov.type, d, state$nF) + 
    state$nH * (state$nF - 1) + 
    switch(state$psi.type,
           EV = state$nH * d,
           EE = d,
           IV = state$nH,
           IE = 1,
           C = 0) + 
    length(state$pi) - 1
  bic = log(nrow(state$zhat)) * df - 2 * state$loglik
  ari = if (is.null(true.classes)) NA else mclust::adjustedRandIndex(rowMaxs(state$zhat), true.classes)
  min.eig = min(apply(state$sigma, 3, function(x) eigen(x, only.values = TRUE)$values))
  
  assign = max.col(state$zhat)
  assign = factor(assign, 1:(nF + nH))
  tab = as.integer(table(assign))
  names(tab) = c(paste0("F", 1:nF), paste0("H", 1:nH))
  
  cat("Fitted Factor-Hybrid Clustering model with",
      nF, "factor and",
      nH, "hybrid clusters.\n")
  
  cat("  BIC (lower is better): ", round(bic, 4), "\n")
  cat("  Log-Likelihood:        ", round(state$loglik, 4), "\n")
  cat("  No. Parameters:        ", df, "\n")
  
  if (!is.null(true.classes))
  cat("  ARI:                   ", round(adjustedRandIndex(assign, true.classes), 4), "\n")
  else if (!is.null(state$true.classes)) 
  cat("  ARI:                   ", round(adjustedRandIndex(assign, state$true.classes), 4), "\n")
  
  cat("\n")
  cat("Component Assignment Counts: \n")
  print(tab)
}

print.FHclust = summary.FHclust