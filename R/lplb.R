#######clean worspace and set file path
#rm(list=ls())
#getwd()

########Load library and functions
library(survival)
library(MASS)

#source('lplb_basicFunctions.R')

lplb <- function(x, ...) UseMethod("lplb")

lplb.default <- function(x, y, control, ...){
  t0 = Sys.time()
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
  Q1 = fit$maxT
  sd = fit$sd
  cat('Q1 = ', Q1, '\n')
  fit$mTstar = bstrp(X, y, control)
  B = control$B
  pvalue = (sum(fit$mTstar>=Q1)+0.5)/(B+1)
  cat('p-value = ', pvalue, '\n')
  t1 = Sys.time()
  runningtime=t1-t0
  fit$w = w
  fit$Q1 = Q1
  fit$B = B
  fit$pvalue = pvalue
  fit$control = control
  fit$call <- match.call()
  #fit$kernel = control$kernel
  #fit$h = control$h
  fit$runningtime = runningtime
  class(fit) <- "lplb"
  
  return(fit)
}

lplb.control = function(h = 0.2, kernel = 'gaussian', B = 200, w0 = seq(0.05, 0.95, 0.025), p1 = 1, pctl = seq(0.2, 0.8, 0.1)) {
  if (!is.numeric(B) || B <= 0) 
    stop("value of 'B' must be > 0")
  if (!is.numeric(h) || h <= 0 || h >= 1) 
    stop("value of 'h' must be in (0, 1)")
  
  if (!is.numeric(p1) || p1 < 1) 
    stop("value of 'p1' must be > or = 1")
  
  #pctl defines the print output of beta(w)
  if (!is.numeric(pctl) || max(pctl) >= 1 || min(pctl)<=0)
    stop("value of 'pctl' must be in (0, 1)")
    
  
  return(list(h = h, B = B, w_est = w0, p1 = p1, pctl = pctl, kernel = kernel))
}

lplb.formula <- function(formula, data=list(...), control = list(...), ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  
  #remove intercept
  x = x[, -1]
  
  y <- model.response(mf)
  
  control = do.call("lplb.control", control)
  
  fit <- lplb.default(x, y, control)
  fit$call <- match.call()
  fit$formula <- formula
  return(fit)
}

print.lplb <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  p1 = ncol(x$beta_w)
  control = x$control
  c_names = as.character(rep(0,p1)) # to save col names
  w_q = quantile(x$w_est, probs = control$pctl, type=3)
  
  opmtrx = matrix(0,p1*2, length(control$pctl))
  for (i in 1:p1){
    opmtrx[(i*2-1), ] = x$beta_w[match(w_q,x$w_est) ,i]
    c_names[(i*2-1)] = colnames(x$beta_w)[i]
    opmtrx[(i*2), ] = x$sd[match(w_q,x$w_est) ,i]
    c_names[(i*2)] = paste0(colnames(x$beta_w)[i],'(sd)')
  }
  
  opmtrx = rbind(w_q, opmtrx)
  rownames(opmtrx) = c('w0',c_names)
  cat("\nCoefficients(w_est quantile):\n")
  print(opmtrx)
  print(x$runningtime)
  cat('p1 =', p1, '; Bootstrap times =', x$B, '\n')
  cat('Kernel type:',x$kernel, '; Bandwidth (h) = ',x$h, '\n')
  cat('Statistic Q1 =', x$Q1, '; p_value =', x$pvalue, '\n')
}

## plot function.
plot.lplb = function(x, scale = c('original', 'transformed'), ...) {
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
