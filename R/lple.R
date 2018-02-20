########Load library and functions
library(survival)
library(MASS)

lple <- function(x, ...) UseMethod("lple")

lple.default <- function(x, y, control, ...){
  X = as.matrix(x)
  ## data needs to be ordered by time
  st = sort(y[, 1], index.return = TRUE)
  idx = st$ix
  y = y[idx, ]
  X = X[idx, ]
  ## transform w into interval (0,1)
  p=ncol(X)-1
  w = X[, p+1]
  X[ ,p+1]=x.cdf(X[ ,p+1])
  X = as.matrix(X)
  
  fit = lple_fit(X, y, control)
  sd = fit$sd
  fit$w = w
  fit$B = B
  fit$control = control
  fit$call <- match.call()
  class(fit) <- "lple"
  
  return(fit)
}

lple.formula <- function(formula, data=list(...), control = list(...), ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  
  #remove intercept
  x = x[, -1]
  
  y <- model.response(mf)
  
  control = do.call("lplb.control", control)
  
  fit <- lple.default(x, y, control)
  fit$call <- match.call()
  fit$formula <- formula
  return(fit)
}

print.lple <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  p1 = ncol(x$beta_w)

  out = cbind(x$w_est, x$beta_w, x$sd)
  colnames(out) = c('w', 'beta(w)', 'sd')
  print(out)
  cat('Kernel type:',x$kernel, '; Bandwidth (h) = ',x$h, '\n')
}

## plot function.
plot.lple = function(x, scale = c('original', 'transformed'), ...) {
  scale = match.arg(scale)
  bw = x$beta_w
  p1 = ncol(bw)
  bw_names = colnames(bw) # keep the col names
  control = x$control
  w_q = quantile(x$w, probs = x$w_est)
  w = switch(scale, original = w_q, transformed = x$w_est)

  if(p1 == 1) {
    par(mfrow = c(1, 1))
    plot(w, bw, xlab = 'w', ylab = 'beta(w)', main = bw_names[1], type = 'l')
  } else if (p1 == 2) {
    par(mfrow = c(1, 2))
  } else if ((p1 <= 4) & (p1 > 2)) {
    par(mfrow = c(2, 2))
  } else if ((p1 <= 6) & (p1 > 4)) {
    par(mfrow = c(2, 3))
  } else {
    p2 = ceiling(sqrt(p1))
    par(mfrow = c(p2, p2))
  }
  if(p1 > 1) {
    for (i in 1:p1) {
      plot(w, bw[, i], xlab = 'w', ylab = 'beta(w)', main = bw_names[i], type = 'l')
    }
  }
}

predict.lple = function(object, newdata, newy = NULL) {
  beta = object$beta_w
  w    = object$w_est

  if(missing(newdata))
    X = object$X
  else
    X = newdata

  p = ncol(X)
  nw = X[, p]
  Z = as.matrix(X[, -p])

  ### approximation beta(w) for w where beta is not estimated
  appxf = function(y, x, xout){ yn=approx(x,y,xout=xout,rule=2)$y}
  bz = apply(beta, 2, appxf, x=w, xout = nw)
 
  xb  = rowSums(Z*bz)
  pred = list(lp = xb, risk = exp(xb))

  if(!is.null(newy)) {
    ## Prediction error for survival data is dfined as -log(Lik) of new data
    status = y[, 2] #y[ ,1] is time; y[ ,2]is status
    if(length(status) != n) 
      stop("Error: new y shall have the same subjects as the new data")

    xb  = predict(fit, X)    #(n*1)
    exb = as.vector(exp(xb))
    pred$pe = -sum(status*(xb - log(cumsum(exb))))
  }
  return(pred)
}
