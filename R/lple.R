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
plot.lple = function(x, ..., scale = c('original', 'transformed')) {
  scale = match.arg(scale)
  bw = x$beta_w
  p1 = ncol(bw)
  bw_names = colnames(bw) # keep the col names
  control = x$control
  w_q = quantile(x$w, probs = x$w_est)
  w = switch(scale, original = w_q, transformed = x$w_est)

  if(p1 == 1) {
    par(mfrow = c(1, 1))
    plot(w, bw, ..., xlab = 'w', ylab = 'beta(w)', main = bw_names[1], type = 'l')
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
      plot(w, bw[, i], ..., xlab = 'w', ylab = 'beta(w)', main = bw_names[i], type = 'l')
    }
  }
}

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

### baseline cumulative hazard function and martingale residuals for LPLE
survfit.lple = function(object, se.fit=TRUE, conf.int = .95) {
  beta = object$beta_w
  w    = object$w_est
  y    = object$y

  ### sort data by time
  idx  = order(y[, 1])
  X    = object$X[idx, ]
  y    = y[idx, ]

  p = ncol(X)
  n = nrow(X)
  nw = X[, p]
  Z = as.matrix(X[, -p])
  nevent = y[, 2]
  time   = y[, 1]
  events = sum(nevent)
  n.risk = n:1

  ### approximation beta(w) for w where beta is not estimated
  bz = apply(beta, 2, .appxf, x=w, xout = nw)
  exb = exp(rowSums(Z*bz))
  rcumsum <- function(x) rev(cumsum(rev(x))) # sum from last to first
  rxb = rcumsum(exb)         #sum over risk set 
  haz = nevent/rxb           #hazard function
  varhaz = nevent/rxb^2      #var for haz
  cumhaz = cumsum(haz)       #Breslow estimate of cumulative hazard
  surv   = exp(-cumhaz)      #survival function S(t)
  residuals = nevent - cumhaz*exb # Martingale residuals
  ### code above has been validated for cumhaz and surv function

  result = list(n = n, events = events, time = time, cumhaz = cumhaz, 
		hazard=haz, surv=surv, n.event = nevent, n.risk = n.risk,
		varhaz=varhaz, residuals = residuals)

  if(se.fit) {
    zval = qnorm(1- (1-conf.int)/2, 0,1)
    varh = cumsum(varhaz)
    std.err = sqrt(varh)
    xx   = ifelse(surv==0,1,result$surv)  #avoid some "log(0)" messages
    tmp1 = ifelse(surv==0, 0, exp(log(xx) + zval*std.err))
    tmp2 = ifelse(surv==0, 0, exp(log(xx) - zval*std.err))
    result = c(result, list(upper = pmin(tmp1, 1), lower = tmp2, 
               std.err=std.err, conf.type='log', conf.int=conf.int))
  }
  class(result) = c('survfit.cox', 'survfit')
  return(result)
  ### see also basehaz()
}

### Residuals of LPLE 
residuals.lple = function(object, type=c("martingale", "deviance")) {
  type = match.arg(type)
  sfit = survfit(object, se.fit = FALSE)
  rr = sfit$residuals
  status = sfit$n.event
  drr = sign(rr)*sqrt(-2*(rr+ifelse(status==0, 0, status*log(status-rr))))
  resid = switch(type, martingale = rr, deviance = drr)
  return(resid)
}
