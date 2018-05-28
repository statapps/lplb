ibs = function(x, ...) UseMethod("ibs")

lple.default <- function(x, ...){
  print("Input must be either coxph or lple object.")
}

ibs.lple= function(object, newdata=NULL, newy = NULL) {
  # When one does not use new data to calculate integrated Brier score, 
  # original data and y will be used 
  y = object$y
  if (is.null(newdata)) {
    newdata = object$X
    newy = y
  } else if(is.null(newy)) 
    stop("To calculate Brier score for newdata, newy cannot be NULL.")

  if(length(newdata[, 1]) != length(newy[, 1])) 
    stop("New data and new y must have same number of observations.")
  
  sf = predict(object, newdata, newy)
  ## Matrix of survival function S(t) rows: subjects, columns: time
  ## where sf$risk = exp(beta(w)*z + g(w))
  S = exp(-sf$risk %*% t(sf$cumhaz))
  time = newy[, 1]
  status = newy[, 2]
  n = length(time)

  Indicator = matrix(0, nrow=n, ncol=n)
  Indicator[lower.tri(Indicator)]=1
  Cens = matrix(status, nrow=n, ncol=n)
  Cens[lower.tri(Cens)]=1

  # fit km curve for censoring
  G_fit = survfit(Surv(time, 1-status)~1)
  G = .appxf(G_fit$surv, x=G_fit$time, xout = time)

  G_M = matrix(G, nrow=n, ncol=n)
  G_M[lower.tri(G_M)]=0
  G_M = G_M+(Indicator %*% diag(G))

  BSM = ((Indicator-S)^2)*Cens/G_M
  BSt = apply(BSM, 2, mean, na.rm = TRUE)
    
  dt = diff(c(0, time))
  return(sum(BSt*dt))
}
