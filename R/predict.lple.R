predict.lple = function(object, newdata, newy = NULL) {
  beta = object$beta_w
  gw   = object$g_w
  w    = object$w_est
  
  sfit = survfit(object, se.fit=FALSE)
  lp   = sfit$lp
  risk = sfit$risk
  if(missing(newdata)) {
    X = object$X
    residuals = sfit$residuals
  }
  else {
    ### if there is new data
    X = newdata
    residuals = NULL
    p = ncol(X)
    n = nrow(X)
    nw = X[, p]
    Z = as.matrix(X[, -p])

    ### approximation beta(w) and g(w) for w where beta and g are not estimated
    bnw = apply(beta, 2, .appxf, x=w, xout = nw)
    gnw = .appxf(gw, x=w, xout = nw)
    lp  = rowSums(Z*bnw + gnw)
    risk = exp(lp)
  }
  result = list(lp = lp, risk = risk, cumhaz = sfit$cumhaz, time = sfit$time)

  if(!is.null(newy)) {
    if(missing(newdata)) stop("Error: newdata cannot missing when there is new y.")
    if(length(newy[, 1]) != n)
      stop("Error: new y shall have the same subjects as the new data.")

    ## sort survival time for new y from smallest to largest
    idx = order(newy[, 1])
    lp  = lp[idx]
    risk = risk[idx]
    newy = newy[idx, ]

    ## Prediction error for survival data is dfined as -log(Lik) of new data
    ## rcumsum is used to calculate at risk observation for ordered time
    time   = newy[, 1]  #time
    status = newy[, 2]  #status
    result$pe = -sum(status*(lp - log(rcumsum(risk))))

    #### prediction error based on martingle residual, may not work as good
    chz = .appxf(sfit$cumhaz, x=sfit$time, xout=time)
    residuals = (status - chz*risk)
    result$pe.mres = sum(residuals^2)
    ## update lp, riks, time and cumhaz, all ordered by time
    result$lp = lp
    result$risk = risk
    result$time = time
    result$cumhaz = chz
  }
  result$residuals = residuals

  return(result)
}
