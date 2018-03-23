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
  ctl = x$control
  cat('Kernel type:',ctl$kernel, '; Bandwidth (h) = ',ctl$h, '\n')
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


##########LPLE fit was checked aganist the following code and get correct results. 
.coxlple = function(X, y, weights) {
  b0 = rep(0, length(X[1, ]))
  ds = 1
  while(ds > 1e-10) {
    sc = .wscore(X, y, b0, weights)
    b0 = b0 + sc 
    ds = sum(abs(sc))
  }
  b0
}

.wscore = function(X, y, beta, weights) {
  status = y[, 2]
  lp = X%*%beta
  exb = exp(lp)
  n = length(status)

  ## Sm is a upper triangular matrix of 1
  Sm = matrix(1, n, n)
  Sm[lower.tri(Sm)] = 0
       
  rcumsum <- function(x) rev(cumsum(rev(x)))
  s0 = rcumsum(exb*weights)
  s1 = Sm%*%((X)*as.vector(exb*weights))
  sc = (X-s1/s0)*as.vector(weights*status)
  V = t(sc)%*%sc
  return(solve(V)%*%t((weights*status)%*%(X-s1/s0)))
}


##############Asymptotic SCB    
.dvn = function(h, kn) {
  ### asymptotic distribution
  A = 5
  ds = 0.001
  s = seq(-A, A, ds)
  
  #normal kernel
  ks = K_func(s, 0, 1, kernel=kn)
  
  ks1=c(0, ks)
  ks2 = c(ks, 0)
  #derivative of kernel
  dks = (ks2-ks1)/ds
  
  # integration of v0 = k(w)^2 and v1 = k'(w)^2
  v0 = mean(ks^2)*2*A
  v1 = mean(dks^2)*2*A
  #cat('v0=', v0, 'v1 = ', v1, '\n')
  dvn = (-2*log(h))^0.5 + (-2*log(h))^-0.5*log(v1/(4*pi*v0))
  return(dvn)
}
#xa = -log(-log(a)/2)
#p=(-2*log(h))^(-0.5) * log((2*pi*log(a^(-0.5)))^-1 ) +  (-2*log(h))^(0.5)

asymSCB = function(h, kernel, conf.int = 0.95){
  dn = .dvn(h, kernel)
  pn = (dn + (log(2)-log(-log(conf.int)))*(-2*log(h))^-0.5)
}
