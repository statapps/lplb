lple_boot = function(object,conf.int = 0.95) {
  X = object$X
  y = object$y
  control = object$control
  B = control$B
  bp = rep(0, B)
  beta_w = object$beta_w
        
  boot_beta = matrix(0, B, length(beta_w))
  cat('\nBootstraping...')
  for(i in 1:B) {
    if(B>30 & (i%%10 == 0)) cat('.')
    idx = sample(1:n, n, replace = TRUE)
    datb = dat[idx, ]
    Xb = cbind(datb$z, datb$w)
    yb = Surv(datb$time, datb$status)
    idx = order(yb[, 1])  # sort by time before call lple_fit
    yb = yb[idx, ]
    Xb = Xb[idx, ]
    fitb = lple_fit(Xb, yb, control)
    beta_b = as.vector(fitb$beta_w)
    sd_b = as.vector(fitb$sd)
    bp[i] = max(abs((beta_b - beta_w)/sd_b))
    boot_beta[i, ] = beta_b
				  }
    boot_se = apply(boot_beta, 2, sd)
    boot_bias = apply(boot_beta, 2, mean) - beta_w
    qbp = quantile(pb, conf.int)
 
    cat('\nDone! Bootstrap scb (p) = ', qbp)
    return(list(boot_beta = boot_beta, bias = boot_bias, se = boot_se, boot_p = bp, qbp = qbp))
}
