
Sk2g = function(X, y, control, bw0, w0, Sm, haz) {
  kernel = control$kernel
  status = y[, 2]
  h = control$h
  w_est = control$w_est
  p1 = control$p1
  X_R = interaction_X_w0(X, p1, w0)
  p = ncol(X) - 1
  p_fai = p + p1 + 1
  n = nrow(X)
  w = X[, ncol(X)]
  xi = as.matrix(bw0)
  
  R = t(X_R)
  H = diag(rep(c(1, h, h), c(p, p1, 1)), p_fai, p_fai)
  H1 = diag(rep(c(1, 1/h, 1/h), c(p, p1, 1)), p_fai, p_fai)
  fai = H  %*% xi
  Tx  = H1 %*% R
  
  kh = K_func(w, w0, h, kernel)
  ezb_fai = as.vector(exp(t(Tx) %*% fai)) * kh
  
  Sf0 = Sm %*% ezb_fai
  Sf1 = Sm %*% (t(Tx) * ezb_fai)
  
  TSf1 = t(Tx) - Sf1/c(Sf0)
  
  TSf1_2 = apply(TSf1, 1, function(x) {return(x %*% t(x))})
  TSf2   = aperm(array(TSf1_2, c(p_fai, p_fai, n)), c(3, 1, 2))

  pi_fai = colSums(status * (TSf2*kh^2), dims = 1)
  B_n    = colSums(kh*TSf1*haz,   dims = 1)
  
  Zt_fai = apply(t(Tx), 1, function(x) {return(x %*% t(x))})
  Z2_fai = array(Zt_fai, c(p_fai, p_fai, n))
  Z2ezb_fai = Z2_fai * array(rep(ezb_fai, each = p_fai * p_fai), 
                          c(p_fai, p_fai, n))
  Sf2 = apply(Z2ezb_fai, c(1, 2), function(x, y) {return(y %*% x)}, Sm)
  Sf1_2 = aperm(array(apply(Sf1, 1, function(x) { return(x %*% t(x))}), 
                    c(p_fai, p_fai, n)), c(3, 1, 2))
  A_n = colSums(status * (kh * (Sf2/c(Sf0) - Sf1_2/c(Sf0)^2)), dims = 1)
  ###I_fai is the same as An in the paper
  return(list(A_n = A_n, pi_fai = pi_fai, B_n = B_n))
}

lple_se = function(X, y, control, betaw, gw){
  h = control$h
  w_est = control$w_est
  m = length(w_est)
  p1 = control$p1
  n = nrow(X)
  p = ncol(X)
  beta = as.matrix(betaw[, 1:p1], n, p1)
  nw = X[, p]
  print(nw)
  nevent = y[, 2]

  Z = as.matrix(X[, -p])

  ### approximation beta(w) and g(w) for w where beta and g are not estimated
  bnw = apply(beta, 2, .appxf, x=w, xout = nw)
  gnw = .appxf(gw, x=w, xout = nw)
  lp  = rowSums(Z*bnw) + gnw
  exb = exp(lp)
  rxb = rcumsum(exb)         #sum over risk set 
  haz = exb*nevent/rxb       #hazard function

  Sm = matrix(1, n, n)
  Sm[lower.tri(Sm)] = 0
  mtrx2 = matrix(0, nrow = p1, ncol = p + p1)
  for (i in 1:p1) mtrx2[i, i] = 1
  
  sd.err = matrix(0, nrow = m, ncol = p1)
  bias = sd.err
  for (i in 1:m) {
    w0 = w_est[i]
    bw0 = betaw[i, ]
    skf = Sk2g(X, y, control, bw0, w0, Sm, haz)
    A_n = skf$A_n
    A1  = solve(A_n)
    pi_fai = skf$pi_fai
    gamma = A1 %*% pi_fai %*% A1 
    gamma11 = mtrx2 %*% gamma %*% t(mtrx2)
    sd.err[i, ] = sqrt(diag(gamma11))
    #print(solve(pi_fai)%*%t(skf$B_fai))
    bias[i, ] = (A1%*%skf$B_n)[1:p1]
  }
  return(list(sd.err = sd.err, bias = bias))
}
