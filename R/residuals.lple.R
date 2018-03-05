### Residuals of LPLE 
residuals.lple = function(object, type=c("martingale", "deviance")) {
  type = match.arg(type)
  sfit = survfit(object, se.fit = FALSE)
  rr = sfit$residuals
  status = sfit$n.event
  drr = sign(rr)*sqrt(-2*(rr+ifelse(status==0, 0, status*log(status-rr))))
  resid = switch(type, martingale = rr, deviance = drr)
  return(resid)
}
