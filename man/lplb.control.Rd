\name{control}
\alias{lplb.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Auxiliary function for lplb fitting}
\description{
Auxiliary function for \code{\link{lplb}} fitting. 
Typically only used internally by 'lplb', but may be used to construct a control argument to either function. 
}
\usage{
#   lplb.control(h, kernel = 'gaussian', B, w0, p1, pctl)
}

\arguments{
  \item{h}{bandwidth of kernel function. The default value is h = 0.2}
  \item{kernel}{kernel funtion types, including "gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine". The default value is 'gaussian'}
  \item{B}{number of bootstrap times. The default value is 200}
  \item{w0}{the estimated points in the interval of (0,1), select arbitrarily. The default value is seq(0.05, 0.95, 0.025)}
  \item{p1}{the number of dependend variables that make interactions with the biomarker w. The default value is 1}
  \item{pctl}{the estimated points that want to be shown in the output. The default value is seq(0.2, 0.8, 0.1)}
}

\details{
Control is used in model fitting of lplb.
}
\value{
This function checks the internal consisitency and returns a list of value as inputed to control model fit of lplb.
}

\author{Siwei Zhang and Bingshu E. Chen}

\seealso{
\code{\link{lplb}}
}
\examples{
## The default control values are: h = 0.2, kernel = 'gaussian', B = 200, w0 = seq(0.05, 0.95, 0.025), p1 = 1, pctl = seq(0.2, 0.8, 0.1) 
##
## To fit the lplb model with some control variables changed, 

w0=seq(0.03,0.97,by=2/100)
ctl = lplb.control(w0=w0, h=0.3, p1=2, B=100)

## then fit the following model

#fit = lplb(x, y, control = ctl)
}
\keyword{ control }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
