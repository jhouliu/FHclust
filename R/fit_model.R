FHclust = function(data, Kf, Kh, cov.type = "VVV", psi.type = "EV", psi.C = 1e-10, init.subsample = 0.5,
                   seed = 12345, mini.em.starts = 100, iters.1 = 100, iters.2 = 10000, AECM.eps = 1e-8, 
                   verbosity = 1, mini.em.min.ev = 1e-4, true.classes = NULL) {

  # Sanity check
  if (!Rfast::is_integer(Kf)) stop("Kf is not an integer!")
  if (!Rfast::is_integer(Kh)) stop("Kh is not an integer!")
  if (Kf < 2) stop("Kf must be at least 2!")
  if (Kh < 1) stop("Kh must be at least 1!")
  psi.type = toupper(psi.type)
  if (length(psi.type) != 1) stop("Invalid psi.type!")
  if (!(psi.type %in% c("EV", "EE", "IV", "IE", "C"))) stop("Unrecognized psi.type!")
  cov.type = toupper(cov.type)
  if (length(cov.type) != 1) stop("Invalid cov.type!")
  cov.types = c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV")
  if (!(cov.type %in% cov.types)) stop("Unrecognized cov.type!")
  if (!is.matrix(data)) tryCatch(data <- as.matrix(data), 
                                 error = function(e) stop("Could not coerce data to matrix."))
  
  if (verbosity > 0) {
    debugmsg = function(...) {message(sprintf(...))}
  } else {
    debugmsg = function(...) {invisible(NULL)}
  }

  # Initialize with mini-EM if so desired
  seed.state = .Random.seed
  set.seed(seed)
  if (iters.1 > 0 & mini.em.starts > 1) {
    debugmsg("%i Mini-EM starts with %i factor and %i hybrid clusters\nwith covariance type %s and error type %s...",
             mini.em.starts, Kf, Kh, cov.type, psi.type)
    states = replicate(mini.em.starts,
                       tryCatch(initialization(data, Kf, Kh, cov.type, psi.type, psi.C = psi.C, subsample = init.subsample),
                                error = function(e) NULL),
                       simplify = FALSE)
    updated.states = lapply(states, function(state) 
      return(tryCatch(FHclust.more.iters(state, iters.1, verbosity > 1, AECM.eps = -1),
                      error = function(e) NULL)))
    state.stats = lapply(updated.states, summary_state, true.classes = NULL)
    state.loglik = sapply(state.stats, `[`, "ll")
    state.min.ev = sapply(state.stats, `[`, "min.eig")
    state.loglik[state.min.ev < mini.em.min.ev] = -Inf
    if (!any(is.finite(state.loglik))) stop("No Mini-EM iterations successful.")
    best.start = which.max(state.loglik)
    state = updated.states[[best.start]]
    rm(states, updated.states)
    debugmsg("Initialization successful.")
  } else {
    debugmsg("Single initialization with %i factor and %i hybrid clusters\nwith covariance type %s and error type %s...",
             Kf, Kh, cov.type, psi.type)
    state = initialization(data, Kf, Kh, cov.type, psi.type, psi.C = psi.C, subsample = init.subsample)
    debugmsg("Initialization successful.")
  }
  .Random.seed = seed.state
  
  # Fit selected model
  if (iters.2 > 0) {
    debugmsg("Beginning %i EM iterations...", iters.2)
    state = FHclust.more.iters(state, iters.2, verbosity > 1, AECM.eps = AECM.eps)
  }

  # Print final message
  if (state$converged) {
    if (iters.1 > 0) {
      debugmsg("Convergence reached at iteration %i + %i.", iters.1, length(state$ll) - iters.1)
    } else {
      debugmsg("Convergence reached at iteration %i.", length(state$ll))
    }
  } else {
    debugmsg("Iteration limit reached without convergence.")
  }
  
  # Attach true class labels if available
  state$true.classes = true.classes
  
  # Set class
  class(state) = "FHclust"
  return(state)
}

FHclust.more.iters = function(state, more.iters = 10000, verbosity = 1, AECM.eps = 1e-2) {
  # if (is.null(state)) return(NULL)
  if (verbosity > 0) {
    debugmsg = function(...) {message(sprintf(...))}
  } else {
    debugmsg = function(...) {invisible(NULL)}
  }
  
  if (more.iters > 0) {
    debugmsg("Running additional iterations...")
    ll = cpp_iteration_fast_progress(state, more.iters, verbosity = verbosity > 1, show_progress = verbosity > 0, AECM_epsilon = AECM.eps)
    state$ll = c(state$ll, ll)
  }
  
  if (state$converged) {
    debugmsg("Convergence reached at iteration %i.", length(state$ll))
  } else {
    debugmsg("Iteration limit reached without convergence.")
  }
  return(state)
}