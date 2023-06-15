weighted.sum = function(list, weights) {
  return(Reduce(`+`, mapply("*", list, weights, SIMPLIFY = FALSE)))
}

mix.epistatic.components = function(mu.list, sigma.list, r) {
  weights = r / sum(r)
  sigma.inv = lapply(sigma.list, solve)
  
  epistatic.sigma = solve(Reduce(`+`, mapply(`*`, weights, sigma.inv, SIMPLIFY = FALSE)))
  epistatic.mu = epistatic.sigma %*% 
    Reduce(`+`, mapply(function(a, b, c) a * b %*% c, weights, sigma.inv, mu.list, SIMPLIFY = FALSE))
  
  return(list(mu = epistatic.mu, sigma = epistatic.sigma))
}

generate.cube.data = function(d, n, lambda.scale = 1) {
  mu.f = t(as.matrix(expand.grid(rep(list(c(-1,1)), d)))) * 1
  parity = apply(mu.f, 2, prod)

  mu.h = mu.f * 2 / 3 / 2
  alpha = apply(mu.h, 2, function(x) {
    z = solve.QP(diag(2^d), rep(0, 2^d), 
                 t(rbind(mu.f, 1, diag(2^d))), 
                 c(x, 1, rep(0, 2^d)), d + 1)
    return(z$solution)
  })
  
  # nullspace = function(x) {
  #   qr = qr(x)
  #   indices = if (qr$rank == 0) 1:ncol(x) else -(1:qr$rank)
  #   return(qr.Q(qr, complete = TRUE)[, indices, drop = FALSE])
  # }
  
  sigma.f = vapply(1:(2^d), function(f) {
    if (parity[f] > 0) {
      0.025 * lambda.scale * diag(d)
    } else {
      mu = mu.f[,f]
      mu = mu / sqrt(sum(mu^2))
      0.018 * lambda.scale * tcrossprod(mu) + 0.002 * lambda.scale * diag(d)
      # rot = cbind(mu, nullspace(mu))
      # rot %*% diag(c(0.02, rep(0.002, d - 1))) %*% t(rot)
    }
  }, matrix(0, d, d))
  
  psi = vapply(1:(2^d), function(f) {
    if (parity[f] > 0) {
      0.001 * lambda.scale * diag(rep_len(c(10, 1), d))
    } else {
      0.001 * lambda.scale * diag(rep_len(c(1, 10), d))
    }
  }, matrix(0, d, d))
  
  sigma.h = vapply(1:(2^d), function(h) {
    weighted.sum(asplit(sigma.f, 3), alpha[,h]^2) + psi[,,h]
  }, matrix(0, d, d))
  
  mu = cbind(mu.f, mu.h)
  
  sigma = array(0, c(d, d, 2 * 2^d))
  sigma[,, 1:(2^d)] = sigma.f
  sigma[,, 2^d + 1:(2^d)] = sigma.h
  sigma_inv = vapply(1:(2*2^d), function(i) solve(sigma[,,i]),
                     matrix(0, d, d))
  
  data.list = lapply(1:(2*2^d), function(g) rmvnorm(n, mu[,g], sigma[,,g]))
  data = do.call(rbind, data.list)
  true.classes = rep(1:(2*2^d), rep(n, (2*2^d)))
  zhat = kronecker(diag(2 * 2^d), rep(1, n))
  
  state = list(data = data, 
               nF = 2^d,
               nH = 2^d,
               zhat = zhat,
               pi = rep(2^-d, 2^d),
               psi.type = "EV",
               cov.type = "VVV",
               stochastic = FALSE,
               hold.z = FALSE,
               convenience = list(Amat = cbind(1L, diag(2^d)),
                                  bvec = c(1L, rep(0, 2^d)),
                                  magic.mat = kronecker(diag(2^d), rep(1, d))),
               true.classes = true.classes,
               mu = mu,
               sigma = sigma,
               sigma_inv = sigma_inv,
               psi = psi,
               alpha = alpha,
               ll = numeric(0),
               loglik = -Inf)
  return(state)
}

generate.cube.data.epistatic = function(d, n, lambda.scale = 1) {
  nF = 2^d
  nH = d * 2^(d-1) 
  mu.f = t(as.matrix(expand.grid(rep(list(c(-1,1)), d)))) * 2
  parity = sign(apply(mu.f, 2, prod))
  sigma.f = vapply(1:(2^d), function(f) {
    mu = mu.f[,f]
    mu = mu / sqrt(sum(mu^2))
    
    if (parity[f] > 0) {
      0.2 * lambda.scale * tcrossprod(mu) + 0.05 * lambda.scale * diag(d)
    } else {
      0.1 * lambda.scale * tcrossprod(mu) + 0.1 * lambda.scale * diag(d)
    }
  }, matrix(0, d, d))
  
  dists = as.matrix(dist(t(mu.f)))
  
  adj.mat = sapply(1:(2^d), function(i) 
    sapply(1:(2^d), function(j) {
      i < j && sum(mu.f[,i] != mu.f[,j]) == 1 
    }))
  
  mixing = t(apply(which(adj.mat, arr.ind = TRUE), 1, function(x) {
    y = numeric(2^d)
    y[x] = 1
    y
  }))
  
  epistatics = lapply(1:nH, function(d) {
    mix.epistatic.components(asplit(mu.f, 2), asplit(sigma.f, 3), mixing[d,])
  })
  
  mu.h = vapply(1:nH, function(d) epistatics[[d]]$mu, numeric(d))
  sigma.h = vapply(1:nH, function(d) epistatics[[d]]$sigma, diag(d))

  mu = cbind(mu.f, mu.h)
  
  sigma = array(0, c(d, d, nF + nH))
  sigma[,, 1:nF] = sigma.f
  sigma[,, nF + 1:nH] = sigma.h
  sigma_inv = vapply(1:(nF + nH), function(i) solve(sigma[,,i]),
                     matrix(0, d, d))
  
  data.list = lapply(1:(nF + nH), function(g) rmvnorm(n, mu[,g], sigma[,,g]))
  
  
  garbo_sigma = diag(3, d)
  garbo_points = rmvnorm(floor(n/2), rep(0, d), garbo_sigma)
  mixing = cbind(0, mixing)
  
  data = do.call(rbind, data.list)
  data = rbind(garbo_points, data)
  
  true.classes = rep(1 + 1:(nF + nH), rep(n, (nF + nH)))
  true.classes = c(rep(1, floor(n/2)), true.classes)
  # zhat = kronecker(diag(nF + nH), rep(1, n))
  zhat = matrix(0, length(true.classes), nF + nH + 1)
  zhat[cbind(1:length(true.classes), true.classes)] = 1
  
  state = list(data = data, 
               nF = ncol(mixing),
               nH = nrow(mixing),
               zhat = zhat,
               pi = rep(1/(nF + nH + 1), (nF + nH + 1)),
               true.classes = true.classes,
               mu = mu,
               sigma = sigma,
               sigma_inv = sigma_inv,
               ll = numeric(0),
               dists = dists,
               alpha = prop.table(mixing, 1),
               loglik = -Inf)
  return(state)
}
