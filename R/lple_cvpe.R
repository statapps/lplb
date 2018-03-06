### prediction error using crossvalidation
lple_cvpe = function(X, y, control, K = 5, Dk = NULL) {
  risk = X[, 1]       # to be used in cIndex

  # sort data by survival time
  idx = order(y[, 1])
  y = y[idx, ]
  X = X[idx, ]
  n = length(X[, 1])

  if(is.null(Dk)) Dk = sample(c(1:K), n, replace = TRUE)
  ctl = control
  pe = 0
  risk = X[, 1]       # to be used in cIndex

  for(i in 1:K) {
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
  # Cross-validation C-index
  cidx = survConcordance(y~risk)$concordance
  return(list(pe=pe, cIndex = cidx))
}
