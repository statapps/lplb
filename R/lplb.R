#######clean worspace and set file path
#rm(list=ls())
#getwd()

########Load library and functions
library(survival)
library(MASS)

#source('lplb_basicFunctions.R')

lplb <- function(x, ...) UseMethod("lplb")

lplb.default <- function(X, y, control, ...)
{
  t0 = Sys.time()
  X = as.matrix(X)
  ## data needs to be ordered by time
  st = sort(y[, 1], index.return = TRUE)
  idx = st$ix
  y = y[idx, ]
  X = X[idx, ]
  ## transform w into interval (0,1)
  p=ncol(X)-1
  X[ ,p+1]=x.cdf(X[ ,p+1])
  X = as.matrix(X)
  
  fit = lple_fit(X, y, control)
  Q1 = fit$maxT
  sd = fit$sd
  cat('Q1 = ', Q1, '\n')
  fit$mTstar = bstrp(X, y, control)
  B = control$B
  pvalue = sum(fit$mTstar>=Q1)/B
  cat('pvalue = ', pvalue, '\n')
  t1 = Sys.time()
  runningtime=t1-t0
  fit$Q1 = Q1
  fit$B = B
  fit$pvalue = pvalue
  fit$control = control
  fit$call <- match.call()
  fit$kernel = control$kernel
  fit$h = control$h
  fit$runningtime = runningtime
  class(fit) <- "lplb"
  
  return(fit)
}

lplb.control = function(h = 0.2, kernel = 'gaussian', B = 200, w_est = seq(0.05, 0.95, 0.025), p1 = 1, pctl = seq(0.2, 0.8, 0.1)) {
  if (!is.numeric(B) || B <= 0) 
    stop("value of 'B' must be > 0")
  if (!is.numeric(h) || h <= 0 || h >= 1) 
    stop("value of 'h' must be in (0, 1)")
  
  if (!is.numeric(p1) || p1 < 1) 
    stop("value of 'p1' must be > or = 1")
  
  #pctl defines the print output of beta(w)
  if (!is.numeric(pctl) || max(pctl) >= 1 || min(pctl)<=0)
    stop("value of 'pctl' must be in (0, 1)")
    
  
  return(list(h = h, B = B, w_est = w_est, p1 = p1, pctl = pctl, kernel = kernel))
}

lplb.formula <- function(formula, data=list(...), control = list(...), ...)
{
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

print.lplb <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  p1 = ncol(x$beta_w)
  control = x$control
  c_names = as.character(rep(0,p1)) # to save col names
  w_q = quantile(x$w_est, probs = control$pctl, type=3)
#   cat("\nCoefficients:\nquantile: 20%    40%    50%    60%    80%\nw_est:",w_q, '\n')
#   for (i in 1:p1){
#     cat(colnames(x$beta_w)[i],': ', x$beta_w[match(w_q,x$w_est) ,i],'\n')
#   }
  ## The code above can not align the texts
  
  opmtrx = matrix(0,p1*2, length(control$pctl))
  for (i in 1:p1){
    opmtrx[(i*2-1), ] = x$beta_w[match(w_q,x$w_est) ,i]
    c_names[(i*2-1)] = colnames(x$beta_w)[i]
    opmtrx[(i*2), ] = x$sd[match(w_q,x$w_est) ,i]
    c_names[(i*2)] = paste0(colnames(x$beta_w)[i],'(sd)')
  }
  
  opmtrx = rbind(w_q, opmtrx)
  rownames(opmtrx) = c('w_est',c_names)
  cat("\nCoefficients(w_est quantile):\n")
  print(opmtrx)
  print(x$runningtime)
  #cat('used time = ',  x$runningtime, '\n')
  cat('p1 =', p1, '; Bootstrap times =', x$B, '\n')
  cat('Kernel type:',x$kernel, '; Bandwidth (h) = ',x$h, '\n')
  cat('Statistic Q1 =', x$Q1, '; p_value =', x$pvalue, '\n')
}

## plot function.
## !!! input must be ordered by w_est!!!!!!!!!
plot.lplb = function(x, ...) {
  w  = x$w_est
  bw = x$beta_w # only the coefs of the interaction terms
  p1 = ncol(bw)
  if(p1 == 1) {
    par(mfrow = c(1, 1))
    plot(w, bw, xlab = 'w', ylab = 'beta(w)', main = 'x1', type = 'l')
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
    bw_names = colnames(bw) # keep the col names
    for (i in 1:p1) {
      plot(w, bw[, i], xlab = 'w', ylab = 'beta(w)', main = bw_names[i], type = 'l')
    }
  }
}

###read data
data1 <- read.table("prostat_with_covar.txt", header=T)
## delete missing data (-9999)
data1=data1[c(-42,-488,-471,-193,-262,-473,-475,-279),]
## order by time. Do this in lplb.default()
# tm = data1$time
# idx    = order(tm)
# data1  = data1[idx,]

## If use lplb.default(), we need to adjust the columns of X to :
# the total number of culumn of X = p+1, p is an integer > = 1.
# column 1 to p of X is treatment and other covariates, such age, etc.
# The last column of X is biomarker w

nameX=names(data1[3:7])
X=cbind(data1[,3],data1[,5:7],data1[,4])
X[, 2] = (X[, 2] - mean(X[, 2]))/sd(X[, 2]) # standarize can help show the shape
# p+1 columns in total: z1 to zp, last column is w
## keep the original colnames of X
names(X)=c(nameX[1],nameX[3:5],nameX[2])

time=data1$time
status=data1$status
y = Surv(time, status) # survival object y

### 2) set controls (parameters)
w_est=seq(0.01,0.99,by=4/100)
# w_est=X[ ,ncol(X)]
# w_est=x.cdf(w_est)
w_est=unique(w_est[order(w_est)]) # w_est must be ordered to plot
h = 0.2 # bandwidth h
p1 = 2
B = 200
control=list(h=h,w_est=w_est,p1=p1, B = B)
# p1: the number of first Zi's with interaction terms with W

mod1 <- lplb(X, y, control)
mod1
plot(mod1)

mod2 = lplb(y~trt+age+weight+pf+ap, data1, B = 3, pctl = seq(0.1, 0.9, 0.1), w_est = w_est,
            kernel = 'ep')
print(mod2)
#plot(mod2)

mod3 = lplb(y~trt+age+weight+pf+ap, data1, control)
print(mod3)
#plot(mod3)

paper1 = lplb(y~z1+z2+z3+z4+w, data=X, w_est=seq(0.01, 0.99, 0.04), p1 = 2, B=200)
print(paper1)
plot(paper1)

paper2 = lplb(y~z1+z2+z3+z4+w, data=X, w_est=seq(0.01, 0.99, 0.01), p1 = 2, B=200)
print(paper2)
plot(paper2)
