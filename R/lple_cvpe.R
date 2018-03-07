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

lple_hSelect = function(X, y, control, step = 0.05, K = 5, faster = FALSE){
  h = seq(0.05, 0.95, step)
  pe = h
  n = length(y[, 1])
  Dk = sample(c(1:K), n, replace = TRUE)

  for (i in 1:length(h)) {
    control$h = h[i]
    pe[i] = lple_cvpe(X, y, control, Dk = Dk, faster=faster)$pe
    cat('h = ', h[i], 'CV Prediction Error =', pe[i], '\n')
  }
  h_opt = h[order(pe)[1]]
  return(list(h_opt = h_opt, h = h, pe = pe))
}
