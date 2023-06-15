forcePD = function(mat, min.eig = 1e-8) {
  if (any(eigen(mat, only.values = TRUE)$values < min.eig)) {
    es = eigen(mat)
    return(es$vectors %*% diag(pmax(es$values, min.eig)) %*% t(es$vectors))
  }
  return(mat)
}

weighted.sum = function(list, weights) {
  return(Reduce(`+`, mapply("*", list, weights, SIMPLIFY = FALSE)))
}

find.outermost.vertices = function(x) {
  # Accepts a matrix with vertices as columns or a list of vectors.
  # Returns indices ordered by decreasing L2 distance from the centre.
  if (is.list(x)) x = do.call(cbind, x)
  centre = rowMeans(x)
  distances = colSums((x - centre)^2)
  return(order(distances, decreasing = TRUE))
}

find.convex.combination = function(ext, relint, min.wt = 0) {
  # Both parameters can be a matrix with vertices as columns or a list of 
  # vectors. Returns a matrix of coefficients that minimize the L2 norm, where 
  # each column is the coefficients corresponding to columns of relint.
  if (is.list(ext)) ext = do.call(cbind, ext)
  if (is.list(relint)) ext = do.call(cbind, relint)
  n.ext = ncol(ext)
  d = nrow(ext)
  quadratic = crossprod(ext)
  if (d < n.ext) {
    # If rank-deficient, we will also minimize x^T x.
    diag(quadratic) = diag(quadratic) + 1
  }
  
  coefs = apply(relint, 2, function(x) {
    solve.QP(quadratic,
             colSums(x * ext),
             cbind(1, diag(n.ext)),
             c(1, rep(min.wt, n.ext)),
             1)$solution
  })
  return(coefs)
}

ncovpar = function(modelname, p, G) {
  switch(
    modelname,
    EII = 1,
    VII = G,
    EEI = p,
    VEI = p + G - 1	,
    EVI = p * G - G + 1,
    VVI = p * G,
    EEE = p * (p + 1) / 2,
    EEV = G * p * (p + 1) / 2 - (G - 1) * p	,
    VEV = G * p * (p + 1) / 2 - (G - 1) * (p - 1),
    VVV = G * p * (p + 1) / 2,
    EVE = p * (p + 1) / 2 + (G - 1) * (p - 1),
    VVE = p * (p + 1) / 2 + (G - 1) * p,
    VEE = p * (p + 1) / 2 + (G - 1),
    EVV = G * p * (p + 1) / 2 - (G - 1),
    -Inf
  )
}
