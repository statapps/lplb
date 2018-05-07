
Sk2g = function(X, y, control, bw0, w0, Haz) {
  kernel = control$kernel
  status = y[, 2]
  h  = control$h
  p1 = control$p1
  XR = interaction_X_w0(X, p1, w0)
  p  = ncol(X) - 1
  p2 = p + p1 + 1 #p2 total vector length for fai
  n  = nrow(X)
  w  = X[, ncol(X)]
  bw = as.matrix(bw0)
  
  H = diag(rep(c(1,   h,   h), c(p, p1, 1)), p2, p2)
  H1= diag(rep(c(1, 1/h, 1/h), c(p, p1, 1)), p2, p2)
  fai= H %*% bw

  U  = XR%*%H1
  
  kh = K_func(w, w0, h, kernel)
  ezb_fai = as.vector(exp(U %*% fai)) * kh
  
  #Sm = matrix(1, n, n)
  #Sm[lower.tri(Sm)] = 0
  #Sf0 = Sm %*% ezb_fai
  #Sf1 = Sm %*% (U * ezb_fai)

  # this code is faster
  Sf0 = rcumsum(ezb_fai)  #reverse cumsum
  ubz = U*ezb_fai
  Sf1 = apply(ubz, 2, rcumsum)

  TSf1 = U - Sf1/c(Sf0)
  
  TSf1_2 = apply(TSf1, 1, function(x) {return(x %*% t(x))})
  TSf2   = aperm(array(TSf1_2, c(p2, p2, n)), c(3, 1, 2))

  Pi_n = colSums(status * (TSf2*kh^2), dims = 1)
  #print(sum(haz))
  B_n  = colSums(kh*TSf1*Haz, dims = 1)
 
  Zt_fai = apply(U, 1, function(x) {return(x %*% t(x))})
  Z2_fai = array(Zt_fai, c(p2, p2, n))
  Z2ezb_fai = Z2_fai * array(rep(ezb_fai, each = p2 * p2), c(p2, p2, n))

  #code with reverse cumsum (rcumsum) runs faster
  #Sf2 = apply(Z2ezb_fai, c(1, 2), function(x, y) {return(y %*% x)}, Sm)
  #print(Sf2[1, , ])
  Sf2 = apply(Z2ezb_fai, c(1, 2), rcumsum)
  #print(Sf2[1, , ])
  Sf1_2 = aperm(array(apply(Sf1, 1, function(x) { return(x %*% t(x))}), 
                    c(p2, p2, n)), c(3, 1, 2))
  A_n = colSums(status * (kh * (Sf2/c(Sf0) - Sf1_2/c(Sf0)^2)), dims = 1)
  return(list(A_n = A_n, Pi_n = Pi_n, B_n = B_n))
}

lple_se = function(X, y, control, betaw, gw){
  h = control$h
  w = control$w_est
  m = length(w)
  p1 = control$p1
  n = nrow(X)
  p = ncol(X)
  beta = as.matrix(betaw[, 1:p1], n, p1)
  nw = X[, p]
  nevent = y[, 2]

  Z = as.matrix(X[, -p])
  #print(Z)

  ### approximation beta(w) and g(w) for w where beta and g are not estimated
  bnw = apply(beta, 2, .appxf, x=w, xout = nw)
  gnw = .appxf(gw, x=w, xout = nw)

  exb = exp(rowSums(Z*bnw) + gnw)
  rxb = rcumsum(exb)     #sum over risk set 
  haz0= nevent/rxb       #hazard function
  haz = exb*haz0
  Haz = exb*cumsum(haz0)

  mtrx2 = matrix(0, nrow = p1, ncol = p + p1)
  for (i in 1:p1) mtrx2[i, i] = 1
 
  sd.err = matrix(0, nrow = m, ncol = p1)
  bias = sd.err
  for (i in 1:m) {
    w0  = w[i]
    bw0 = betaw[i, ]
    skf = Sk2g(X, y, control, bw0, w0, Haz)
    A_n = skf$A_n
    A1  = solve(A_n)
    Pi_n= skf$Pi_n
    gamma  = A1 %*% Pi_n %*% A1 
    gamma11= mtrx2 %*% gamma %*% t(mtrx2)
    sd.err[i, ]= sqrt(diag(gamma11))
    bias[i, ]  = (A1%*%skf$B_n)[1:p1]
  }
  return(list(sd.err = sd.err, bias = bias, haz = haz0))
}
