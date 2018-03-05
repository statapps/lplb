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
