lple_boot = function(object, method = c("residual", "simple"), 
  conf.int = 0.95, parallel = TRUE) {
  method = match.arg(method)

  X = object$X
  y = object$y
  control = object$control
  B = control$B
  status = y[, 2]
  n = length(status)
  pb = rep(0, B)
  beta_w = object$beta_w
  boot.beta = matrix(0, B, length(beta_w))
 
  if(method == "residual") {  
    prd = predict(object)
    mrsd = prd$residuals    # Martingle residuals
    exb  = prd$risk         # sf$risk = exp(b(w)*z)
  }

  idx = matrix(sample(1:n, n * B, replace = TRUE), B, n)
  cat("\nBootstraping...\n")

  fn = function(r) {
    id = idx[r, ]
    
    if(method == "residual") {
      status_sample = status[id]
      resid_sample = status_sample - mrsd[id]
      lambda_star = resid_sample/exb
      
      ord = order(lambda_star)
      Xstar = X[ord, ]
      time_star = 1:n
      status_star = status_sample[ord]
      ystar = Surv(time_star, status_star)
    } else if(method == "simple") {
      Xb = X[id, ]
      yb = y[id, ]
      ord = order(yb[, 1])
      Xstar = Xb[ord, ]
      ystar= yb[ord, ]
    }
    
    fitb = lplb::lple_fit(Xstar, ystar, control)
    beta_b = as.vector(fitb$beta_w)
    sd_b = as.vector(fitb$sd)
    return(list(beta_b, sd_b))
  }
  
  if(parallel) {
    n.cores = detectCores() - 1
    cl = makeCluster(n.cores)
    clusterEvalQ(cl, library(survival))
    cat("Note: Parallel with", n.cores, "cores.\n")
  
    #X
    #status
    #control
    #idx
    boots = parLapply(cl, seq_len(B), fn)
    stopCluster(cl)
  } else {
    boots = lapply(seq_len(B), fn)
  }
  
  for (r in seq_len(B)) {
    beta_b = boots[[r]][[1]]
    boot.beta[r, ] = beta_b
    pb[r] = max(abs(beta_b - beta_w)/boots[[r]][[2]])
  }
  boot.se = apply(boot.beta, 1, sd)
  #qbp = quantile(pb, conf.int)
  #cat("\nDone! Bootstrap scb (p) = ", qbp, "\n")
  structure(list(boot.beta = boot.beta, boot.se = boot.se, 
	#boot.p = pb, qbp = qbp,
	fit = object, method = method), class = "lple_boot")
}

plot.lple_boot = function(x, ...) {
  fit = x$fit
  control = fit$control
  w=control$w_est
  B = length(x$boot.beta[, 1])
  
  plot(w, fit$beta_w, type = 'n', ylim = c(min(x$boot.beta), max(x$boot.beta)))
  for(i in 1:B)
    lines(w, x$boot.beta[i, ], lty = 2, col = 'grey')
  lines(w, fit$beta_w, type = 'l', lwd = 2, col = 'blue')
}
