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
