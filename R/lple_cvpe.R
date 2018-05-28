### prediction error using crossvalidation
lple_cvpe = function(X, y, control, K = 10, Dk = NULL) {
  risk = X[, 1]       # to be used in cIndex
  ibs = 0             # to be used in integrated Brier Score
  pe  = 0             # to be used in likelihood goodness of fit predErr
  # sort data by survival time
  idx = order(y[, 1])
  y = y[idx, ]
  X = X[idx, ]
  n = length(X[, 1])

  if(is.null(Dk)) Dk = sample(c(1:K), n, replace = TRUE)
  else K = max(Dk)

  ctl = control
  pe = 0
  risk = X[, 1]       # to be used in cIndex

  for(i in 1:K) {
    Xi = X[Dk != i, ]
    yi = y[Dk != i, ]
	        
    X.test = X[Dk == i, ]
    y.test = y[Dk == i, ]
   
    fit = lple_fit(Xi, yi, ctl, se.fit=FALSE)
   
    prd = predict(fit, X.test, y.test)
    ibs = ibs + ibs(fit, X.test, y.test)
    # risk score and prediction error
    risk[Dk == i] = prd$risk
    pe = pe + prd$pe
  }
  ibs = ibs/K
  # Cross-validation C-index
  cidx = survConcordance(y~risk)$concordance
  return(list(cIndex = cidx, ibs = ibs, pe = pe))
}

lple_hSel = function(X, y, control, method = c("ibs", "pe"), m = 28, K = 10, parallel = TRUE){
  #hx = seq(0.025, 0.5, step)
  X = X
  y = y
  method = match.arg(method)
  control = control

  n = length(y[, 1])
  alpha = seq(2, 5, length.out=m)
  hx = n^(-1/alpha)
  pe = hx
  sqm = seq_len(length(hx))
  
  Dk = sample(c(1:K), n, replace = TRUE)

  fn = function(r) {
    control$h = hx[r]
    lpe = lple_cvpe(X, y, control, Dk = Dk)
    pe = switch(method, ibs = lpe$ibs, pe = lpe$pe)
  }

  if(parallel) {
    ##Number of cores
    n.cores = detectCores() - 1
    cl = makeCluster(n.cores)
    cat('Note: Parallel with ', n.cores, 'cores.\n')
    #clusterExport(cl, c("X", "y", "control"))
    #clusterEvalQ(cl, library(lplb))
    pe = parSapply(cl, sqm, fn)
    stopCluster(cl)
  } else 
    pe = sapply(sqm, fn)
  
  h_opt = hx[order(pe)[1]]
  a = -log(n)/log(h_opt)
  cat('Optimal h = power(n, -1/', a, ') = ', h_opt, '\n')
  return(list(h_opt = h_opt, h = hx, pe = pe))
}
