### prediction error using crossvalidation
lple_cvpe = function(X, y, control, K = 5, Dk = NULL, faster = FALSE) {
  risk = X[, 1]       # to be used in cIndex

  # sort data by survival time
  idx = order(y[, 1])
  y = y[idx, ]
  X = X[idx, ]
  n = length(X[, 1])

  if(is.null(Dk)) Dk = sample(c(1:K), n, replace = TRUE)
  else K = max(Dk)

  if (faster) m = 1
  else m = K 

  ctl = control
  pe = 0
  risk = X[, 1]       # to be used in cIndex

  for(i in 1:m) {
    Xi = X[Dk != i, ]
    yi = y[Dk != i, ]
	        
    X.test = X[Dk == i, ]
    y.test = y[Dk == i, ]
   
    fit = lple_fit(Xi, yi, ctl)
   
    prd = predict(fit, X.test, y.test)
    # risk score and prediction error
    risk[Dk == i] = prd$risk
    pe = pe + prd$pe
  }
  if(faster) pe = pe * K
  # Cross-validation C-index
  cidx = survConcordance(y~risk)$concordance
  return(list(pe=pe, cIndex = cidx))
}

lple_hSel = function(X, y, control, step = 0.025, K = 5, parallel = c("yes", "no")){
  parallel = match.arg(parallel)
  hx = seq(0.05, 0.95, step)
  pe = hx
  sqm = seq_len(length(hx))
  n = length(y[, 1])
  Dk = sample(c(1:K), n, replace = TRUE)

  fn = function(r) {
    control$h = hx[r]
    pe = lple_cvpe(X, y, control, Dk = Dk)$pe
  }

  if(parallel == "no") {
    pe = sapply(sqm, fn)
  } else {
    ##Number of cores
    n.cores = detectCores() - 1
    cl = makeCluster(n.cores)
    cat('Note: Parallel with ', n.cores, 'cores.\n')
    X
    y
    control
    hx
    #clusterExport(cl, c("X", "y", "control"))
    #clusterEvalQ(cl, library(lplb))
    pe = parSapply(cl, sqm, fn)
    stopCluster(cl)
  }
  h_opt = h[order(pe)[1]]

  a = -log(n)/log(h_opt)
  if ((2 < a) & (a < 5))
    cat('h = power(n, -1/', a, ') = ', h_opt, '\n')
  else {
    cx = h_opt/n^(-1/3.5)
    cat('h = ', cx, '* power(n, -1/3.5) =', h_opt, '\n')
  }
  return(list(h_opt = h_opt, h = h, pe = pe))
}
