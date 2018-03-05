predict.lple = function(object, newdata, newy = NULL) {
  beta = object$beta_w
  w    = object$w_est
  sfit = survfit(object, se.fit=FALSE)

  if(missing(newdata)) {
    X = object$X
    residuals = sfit$residuals
  }
  else {
    X = newdata
    residuals = NULL
  }
  p = ncol(X)
  n = nrow(X)
  nw = X[, p]
  Z = as.matrix(X[, -p])

  ### approximation beta(w) for w where beta is not estimated
  bz = apply(beta, 2, .appxf, x=w, xout = nw)
  xb  = rowSums(Z*bz)
  exb = exp(xb)
  result = list(lp = xb, risk = exb)

  if(!is.null(newy)) {
    if(length(newy[, 1]) != n)
      stop("Error: new y shall have the same subjects as the new data")

    ## sort survival time for new y
    idx = order(newy[, 1], decreasing=TRUE)
    xb  = xb[idx]
    exb = exb[idx]
    newy = newy[idx, ]

    ## Prediction error for survival data is dfined as -log(Lik) of new data
    time   = newy[, 1]
    status = newy[, 2]             #newy[ ,1] is time; newy[ ,2]is status
    result$pe = -sum(status*(xb - log(cumsum(exb))))

    #### prediction error based on martingle residual, may not work as good
    chz = .appxf(sfit$cumhaz, x=sfit$time, xout=time)
    residuals = (status - chz*exb)
    result$pe.mres = sum(residuals^2)
  }
  result$residuals = residuals
  return(result)
}
