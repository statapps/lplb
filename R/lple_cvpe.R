### prediction error using crossvalidation
lple_cvpe = function(X, y, control, K = 10, Dk = NULL, faster = FALSE) {
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

lple_hSel = function(X, y, control, m = 28, K = 10, parallel = TRUE){
  #hx = seq(0.025, 0.5, step)
  X = X
  y = y
  control = control
  n = length(y[, 1])
  alpha = seq(2, 5, length.out=m)
  hx = n^(-1/alpha)
  pe = hx
  sqm = seq_len(length(hx))
  
  Dk = sample(c(1:K), n, replace = TRUE)

  fn = function(r) {
    control$h = hx[r]
    pe = lple_cvpe(X, y, control, Dk = Dk)$pe
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
