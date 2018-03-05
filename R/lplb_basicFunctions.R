### 01. basic functions
## 01_1) transform w into interval (0,1)
x.cdf = function(x) {
  n = length(x)
  p = rep(0, n)
  for (i in 1:n)
    p[i] = sum(x<=x[i])
  p = p/(n+1)
  return(p)
}

## 01_2) Kernel function
K_func<-function(w, u, h, kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")) {
  kernel = match.arg(kernel)
  x = w-u
  ax = abs(x)
  esp = 1e-40

  kh = switch(kernel, gaussian = ifelse(ax < 5*h, dnorm(x, sd = h/2), esp), # I would set the default for guassian
         rectangular = ifelse(ax < h, 0.5/h, esp), 
         triangular = ifelse(ax < h, (1 - ax/h)/h, esp),
         epanechnikov = ifelse(ax < h, 3/4 * (1 - (ax/h)^2)/h, esp), ## This was the kernel that we used before
         biweight = ifelse(ax < h, 15/16 * (1 - (ax/h)^2)^2/h, esp),
         cosine = ifelse(ax < h, (1 + cos(pi * x/h))/(2*h), esp),
         optcosine = ifelse(ax < h, pi/4 * cos(pi * x/(2*h))/h, esp)
         )
  # in previous version, only epanechnikov kernel can be used.
  # kh = (3*(1-((w-u)/h)^2)/4*(abs((w-u)/h)<1)+(abs((w-u)/h)>=1)*1e-40)/h
  return(kh)
}

## 01_3) transform matrix X into new matrix with p1 interaction (with w) terms
interaction_X=function(X, p1){
  p=ncol(X)-1
  n=nrow(X)
  X2 = matrix(0, n, p+p1+1)
  w = X[, p+1]
  X1 =X[, 1:p1]
  # Assign value to X2
  X2[, 1:p] = X[, 1:p]
  X2[, (p+1):(p+p1)] = diag(w)%*%X1
  X2[, p+p1+1] = w
  return(X2)
}
## 01_4) transform matrix X into new matrix with p1 interaction with (w-w0) terms
interaction_X_w0=function(X, p1, w0){
  p=ncol(X)-1
  n=nrow(X)
  ## build interaction terms with (w-w0) in X
  X2 = matrix(0, n, p+p1+1)
  w = X[ ,p+1]
  X1 =X[ ,1:p1]
  X2[ ,1:p] = X[ ,1:p]
  X2[ ,(p+1):(p+p1)] = diag(w-w0)%*%X1 # (n*p1)
  X2[, p+p1+1] = w-w0
  return(X2)
}

### Approximate function
.appxf = function(y, x, xout){ approx(x,y,xout=xout,rule=2)$y }

### 02. Calculate S(k)_theta and S(k)_fai
# Since St0 to St2 do not change w0, so do I_theta and pi_theta, we can separate Sk2 into two parts: Sk2t and Sk2f
# including the calculation of I and pi
## 02_1 S(k)_theta and I_theta, pi_theta
Sk2t = function(X, y, theta) {
  ## NOTE!!!!: input X should be ordered by time already!!!!!!!!!!!!!
  ## setup parameters
  status = y[, 2] #y[ ,1] is time; y[ ,2]is status
  p_th = ncol(X)
  n = nrow(X)
  Gx = t(X)
  ## Sm is a upper triangular matrix of 1
  Sm = matrix(1, n, n)
  Sm[lower.tri(Sm)] = 0

  ezb_th = as.vector(exp(t(Gx)%*%theta)) # (n*1)
  St0 = Sm%*%ezb_th # (n*1)
  St1 = Sm%*%(t(Gx)*ezb_th) # (n*p_th)

  GSt1 = t(Gx)-St1/c(St0) # (n*p_th)
  # Gi - S1/S0, we also need this to calculate pi_cov
  GSt1_2= apply(GSt1, 1, function(x){return(x%*%t(x))}) # ((p_th*p_th)*n)
  GSt2 = aperm(array(GSt1_2, c(p_th, p_th, n)), c(3, 1, 2))
  # aperm changes a p1*p2*p3 array to a p3*p1*p2 array with c(3, 1, 2)
  # GSt2 is a n*p*p array with ith pxp matrix = (Gi-S1i/S0i)%*%t(Gi-S1i/S0i)
  pi_theta = colSums(status*GSt2, dims = 1)

  Zt_th = apply(t(Gx), 1, function(x){return(x%*%t(x))})
  Z2_th = array(Zt_th, c(p_th, p_th, n))
  # multiply each Z2(p, p, i) with ezb[i], by change ezb to a p*p*n array with each of i th pxp matrix = ezb[i]
  Z2ezb_th = Z2_th * array(rep(ezb_th, each = p_th*p_th), c(p_th, p_th, n))
  ## calculate S2, a n*p*p array
  St2 = apply(Z2ezb_th, c(1, 2), function(x, y){return(y%*%x)}, Sm)
  # aperm changes a p1*p2*p3 array to a p3*p1*p2 array with c(3, 1, 2)
  St1_2 = aperm(array(apply(St1,1,function(x){return(x%*%t(x))}),c(p_th,p_th,n)),c(3,1,2))
  I_theta = colSums(status*(St2/c(St0)-St1_2/c(St0)^2), dims = 1)
  # need GSt1= Gi-S1i/S0i for pi_cov
  return(list(I_theta = I_theta, pi_theta=pi_theta, GSt1=GSt1))
}

## 02_2 S(k)_fai and I_fai, pi_fai, pi_cov
Sk2f = function(X, y, control, bw0, w0, GSt1) {
  kernel = control$kernel
  ## NOTE!!!!: input X should be ordered by time already!!!!!!!!!!!!!
  ## setup parameters
  status = y[, 2] #y[ ,1] is time; y[ ,2]is status
  h = control$h
  w_est = control$w_est
  p1 = control$p1

  X_R = interaction_X_w0(X,p1,w0)
  p = ncol(X)-1
  p_th = ncol(X)
  p_fai = p+p1+1
  n = nrow(X)
  w = X[ ,ncol(X)]

  xi = as.matrix(bw0) # class: numeric (a row/col of a matrix)
  R = t(X_R) # (p_fai * n); depend on w0
  H = diag(rep(c(1,h,h),c(p,p1,1)),p_fai,p_fai) # (p_fai * p_fai)
  fai = H%*%xi # (p_fai * 1)
  Tx = solve(H)%*%R # (p_fai * n); depend on w0
  # Tx = (T1, T2, ..., Tn), each Ti is a column vector of (p_fai*1)

  ## Sm is a upper triangular matrix of 1
  Sm = matrix(1, n, n)
  Sm[lower.tri(Sm)] = 0
  ## setup Kernel weight vector: kh
  kh=K_func(w, w0, h, kernel)

  ezb_fai = as.vector(exp(t(Tx)%*%fai)) * kh # (n*1)
  Sf0 = Sm%*%ezb_fai  # vector: (n*1)
  Sf1 = Sm%*%(t(Tx)*ezb_fai) # matrix: (n*p_fai)

  TSf1 = t(Tx)-Sf1/c(Sf0) # (n*p_fai)
  TSf1_2= apply(TSf1, 1, function(x){return(x%*%t(x))}) # ((p_fai*p_fai)*n)
  TSf2 = aperm(array(TSf1_2, c(p_fai, p_fai, n)), c(3, 1, 2))
  # TSf2 is a n*p*p array with ith pxp matrix = (Ti-S1i/S0i)%*%t(Ti-S1i/S0i)
  pi_fai = colSums(status*(TSf2*kh^2), dims = 1)

  Zt_fai = apply(t(Tx), 1, function(x){return(x%*%t(x))})
  Z2_fai = array(Zt_fai, c(p_fai, p_fai, n))
  Z2ezb_fai = Z2_fai * array(rep(ezb_fai, each = p_fai*p_fai), c(p_fai, p_fai, n))
  ## calculate Sf2, a n*p*p array
  Sf2 = apply(Z2ezb_fai, c(1, 2), function(x, y){return(y%*%x)}, Sm)
  Sf1_2 = aperm(array(apply(Sf1,1,function(x){return(x%*%t(x))}),c(p_fai,p_fai,n)),c(3,1,2))
  I_fai = colSums(status*(kh*(Sf2/c(Sf0)-Sf1_2/c(Sf0)^2)), dims = 1)
  ## calculate pi_cov
  # Note: TSf1 = t(Tx)-Sf1/c(Sf0) # (n*p_fai)
  # Note: GSt1 = t(Gx)-St1/c(St0) # (n*p_th)
  pi_cov = t(TSf1)%*%diag(status*kh)%*%GSt1 #p_fai*p_th
  return(list(I_fai = I_fai, pi_fai=pi_fai, pi_cov=pi_cov))
}

### 03. Calculate maximum of Q1 (for 1:m), m is the length of w_est
maxTest = function(X,y,control,theta, betaw){
  ## NOTE!!!!: input X, status should be ordered by time already!!!!!!!!!!!!!
  ## set parameters
  h = control$h
  w_est = control$w_est
  m = length(w_est)
  p1 = control$p1
  n=nrow(X)
  p=ncol(X)-1
  w=X[ ,p+1]

  ## build matrix of beta(w)-beta0, col number: p1
  beta_hat = matrix(rep(theta[1:p1],time=m),nrow=m,ncol=p1,byrow=T)
  betaw_hat = betaw[, 1:p1]
  diff_est = betaw_hat-beta_hat

  ## build var(beta(w)-beta)
  ## I_theta, pi_theta, sigma do not change with w0, so do St0, St1, St2, which should be outside the loop m
  skt=Sk2t(X, y, theta)
  I_theta = skt$I_theta
  pi_theta = skt$pi_theta
  GSt1 = skt$GSt1 # for calculating pi_cov in Sk2f()
  sigma=solve(I_theta)%*%pi_theta%*%solve(I_theta)
  ## constuct return matrixs: sigma11
  mtrx1=matrix(0,nrow=p1,ncol=p+1)
  for (i in 1:p1) mtrx1[i,i]=1
  sigma11=mtrx1%*%sigma%*%t(mtrx1)

  mtrx2=matrix(0,nrow=p1,ncol=p+p1+1)
  for (i in 1:p1) mtrx2[i,i]=1
  # for later use in loop to constract gamma11 and omiga11

  ## initial
  Q1 = rep(0,time=m)
  sd.err = matrix(0,nrow=m,ncol=p1) # to save var, then calculate std.err
  ## loop on w_est
  for (i in 1:m){
    w0  = w_est[i]
    bw0 = betaw[i, ]
    skf=Sk2f(X, y, control, bw0, w0, GSt1)
    I_fai = skf$I_fai
    pi_fai = skf$pi_fai
    pi_cov = skf$pi_cov
    gamma=solve(I_fai)%*%pi_fai%*%solve(I_fai)
    omiga=solve(I_fai)%*%pi_cov%*%solve(I_theta)
    ## constuct return matrixs: gamma11, omiga11
    gamma11=mtrx2%*%gamma%*%t(mtrx2)
    omiga11=mtrx2%*%omiga%*%t(mtrx1)
    
    ## save diag(gamma11) to calculate sd.err of beta_w
    sd.err[i, ] = sqrt(diag(gamma11))

    ## calculate Q1 at w0=w_est[i]
    var_hat = sigma11 + gamma11 - 2*omiga11
    diff_hat = diff_est[i, ] # class: numeric
    Q1[i] = as.numeric( t(diff_hat) %*% solve(var_hat) %*% as.matrix(diff_hat))
  }
  maxQ1 = max(abs(Q1))
  return(list(maxQ1=maxQ1, sd.err=sd.err))
}

### 04. estimate beta(w) (Under H1) and theta (Under H0)
lple_fit = function(X, y, control) {
  h = control$h
  kernel = control$kernel
  w_est = control$w_est
  x_names = colnames(X)
  p1 = control$p1    # number of interaction terms between Zi and w
  p = length(X[1, ])-1 # number of Zi's
  n = length(X[ ,1]) # number of observations
  m = length(w_est) # number of estimate points; choose from intervel (0,1) arbitrarily
  X = as.matrix(X)  # coxph() need X to be class matrix

  ## X_fai is the data of the interaction model (H1), estimator is beta(w)
  X_fai=interaction_X(X,p1)
  ## X is the data of the no-interaction model (H0), estimator is theta
  fitH0=coxph(y ~ X)
  theta=fitH0$coef

  ## estimate beta(w)
  ## set matrix "betaw" to save the coef of Z1 to Zp, 
  ## interaction terms (zi with (w-w0)) and (w-w0) in every row, 
  ## each row for a different w0. dim: (m* (p+p1+1))
  betaw = matrix(0, nrow=m, ncol=p+p1+1)
  ## set matrix "beta_w" to save only the Z1 to Zp1, 
  ## by which we can plot the relationship of beta(w) vs. w_est. dim: (m*p1)
  beta_w = matrix(0, m, p1)

  ## Check if the data is sorted by time
  j = sample(1:n, 2)
  tm = y[j, 1]
  if((j[1]-j[2])*(tm[1]-tm[2])<0) 
    stop("Survival time shall be sorted! Check your program or use lple() directly")
  for (i in 1:m) {
    wg = K_func(X_fai[ ,p+p1+1], w_est[i],h, kernel)
    fit = coxph(y ~ X_fai, subset= (wg>0), weights=wg)
    betaw[i, ] = fit$coef
  }
  beta_w[, 1:p1] = betaw[ ,(1:p1)]+betaw[ ,(p+1):(p+p1)] * w_est
  betaw[, 1:p1]= beta_w

  ## return value
  maxreturn = maxTest(X, y, control, theta, betaw)
  maxT = maxreturn$maxQ1
  sd = maxreturn$sd.err
  colnames(beta_w) = x_names[1:p1]
  fit = list(w_est = w_est, beta_w = beta_w, betaw = betaw, maxT = maxT, 
	     sd=sd, h = h, X = X, y = y, kernel = kernel)
  class(fit)= "lple"
  fit$call = match.call()
  return(fit)
}

### 05. bootstrap
bstrp = function(X, y, control){
  X = as.matrix(X)
  kernel = control$kernel
  h  = control$h
  p1 = control$p1
  B  = control$B
  ## Fit model under H0
  fitH0=coxph(y ~ X)
  exb  = exp(X%*%fitH0$coef) #to be used for residual bootstrap

  X_fai = interaction_X(X, p1)
  w = X[ ,ncol(X)]
  n=nrow(X)
  status = y[, 2] #y[ ,1] is time; y[ ,2]is status, to be the initial value

  ## build Martingale residuals
  resid_1 = rep(0, n)
  for (i in 1:n) {
    wg = K_func(w, w[i], h, kernel) #do not use w_est here, but use w[i], length(w) = n
    fit = coxph(y ~ X_fai, subset= (wg>0), weights=wg)
    resid_1[i]<-residuals(fit, "martingale")[i]
  }

  # Bootstraping
  mTstar = rep(0, B)
  i = 1
  cat('Bootstraping')
  Bn = floor(B/10) + 1     
  while (i <= B) {
    sample_index  = sample(n, size=n,replace=T)
    status_sample = status[sample_index]
    resid_sample  = status_sample - resid_1[sample_index]
    
    # create new survival time T_star using lambda_star
    lambda_star = resid_sample/exb

    # sort dataset by lambda_star
    time_order  = order(lambda_star)
    status_star = status_sample[time_order]
    X_star      = X[time_order, ]
    time_star   = 1:n
    y_star      = Surv(time_star, status_star)

    g = try(lple_fit(X_star,y_star,control))
    # use another set of bootstrap if coxph does not converge within 100 iter
    if(class(g) == "try-error") next
    # next halts the processing of the current iteration and advances the looping index.
    mTstar[i] = g$maxT
    #cat('i = ', i, mTstar[i], '\n') # for inspection
    #cat(' ', i, ',',sep = "") # for simplified inspection
    if((i %% Bn)==0) cat('.')
    i = i + 1
  }
  cat('DONE!\n')
  return(mTstar)
}
