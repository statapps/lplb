lple_bootp = function (object, conf.int = 0.95) {
  X = object$X
  y = object$y
  control = object$control
  B = control$B
  pb = rep(0, B)
  beta_w = object$beta_w
  boot.beta = matrix(0, B, length(beta_w))
  cat('\nBootstraping...\n')

  ##Number of cores
  n.cores = detectCores() - 1
  cl = parallel::makeCluster(n.cores)
  
  cat('number of cores = ', n.cores, '\n')
  #clusterExport(cl, c("X", "y"))
  #parallel::clusterEvalQ(cl, library(lplb))
  idx = matrix(sample(1:n, n*B, replace = TRUE), B, n)
  X
  y  
  control
  fn = function(r) {
    id = idx[r, ]
    Xb = X[id, ]
    yb = y[id, ]
    ord = order(yb[, 1])
    yb = yb[ord, ]
    Xb = Xb[ord, ]
    fitb = lple_fit(Xb, yb, control)
    beta_b = as.vector(fitb$beta_w)
    sd_b = as.vector(fitb$sd)
    return(list(beta_b, sd_b))
  }
  ## resampling using parallel computing
  res = parallel::parLapply(cl, seq_len(B), fn)
  parallel::stopCluster(cl)

  ## Clean up the data
  for (r in seq_len(B)) {
    beta_b = res[[r]][[1]]
    boot.beta[r, ] = beta_b
    pb[r] = max(abs(beta_b - beta_w)/res[[r]][[2]])
  }
  boot.se = apply(boot.beta, 2, sd)
  qbp = quantile(pb, conf.int)
  cat("\nDone! Bootstrap scb (p) = ", qbp, '\n')
  return(list(boot.beta = boot.beta, boot.se = boot.se, boot.p = pb, 
        qbp = qbp))
}
