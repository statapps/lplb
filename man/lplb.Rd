\name{lplb}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{lplb}
\alias{lple_fit}
\alias{bstrp}
\alias{lplb.default}
\alias{lplb.formula}

\title{
Local partial likelihood bootstrap (LPLB) method to fit biomarker Models
}
\description{
\{lplb\} is a R package for local partial likelihood estimation (LPLE) (Fan et al., 2006) of the coefficients of covariates with interactions of the biomarker W, and hypothesis test of whether the relationships between covariates and W are significant, by using bootstrap method.
}
\usage{
\method{lplb}{default}(x, y, control, ...)
\method{lplb}{formula}(formula, data=list(...), control = list(...), ...)

# use
#          lplb(y ~ X1+X2+...+Xp+w, data=data, control)
#
# to fit a model with interactions between biomarker (w) with the first p1
# terms of dependent variables. 
# p1 is included in 'control'. p1<p. See 'lplb.control' for details
#
# use
#         lplb(x, y, control)
#
# to fit a model without the formula
#
# Biomarker w should be the 'LAST' dependend variable
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class "formula"(or one that can be coerced to that class): a symbolic description of the model to be fitted.
  The details of model specification are given under 'Details'.}

  \item{data}{an optional data frame, list or environment (or object coercible by 'as.data.frame' to a data frame) containing the variables in the model.
If not found in data, the variables are taken from environment(formula).}

  \item{x, y}{For 'lplb.default', x is a design matrix of dimension n * (p+1) and y is a vector of observations of length n for a "Surv" object for "coxph".}

  \item{control}{a list of parameters for controlling the fitting process. See 'lplb.control' for details}

  \item{\dots}{additional arguments to be passed to the low level regression fitting functions (see below).} 
%- do we need this??????????
}
\details{
Here 'w' is a Biomarker variable. This variable is required and shall be the last dependent variable in the formula.

'x.cdf' is a function that maps biomarker values to interval (0, 1) using its empirical cumulative distribution function.
}
\value{
%%  ~Describe the value returned
lplb returns an object of class inheriting from "lplb" which inherits from the class 'coxph'. See later in this section.

The function "print" (i.e., "print.lplb") can be used to obtain or print a summary of the results.

An object of class "lplb" is a list containing at least the following components:

  \item{beta_w}{a matrix of m * p1, the estimated coefficients at each of the m estimated points, for the first p1 dependent variables with interactions of the biomarker w}
  \item{Q1}{the test statistic of the data}
  \item{mTstar}{a vector of the test statistics from B times' bootstrap}
  \item{pvalue}{the p-value of the hypothesis that beta_w is a constant}
}
\references{
Zhang, S., Jiang, W. and Chen, B. E. (2016). Estimate and test of multivariate covariates 
and biomarker interactions for survival data based on local partial likelihood. Manuscript in preparation. 
}
\author{
Siwei Zhang and Bingshu E. Chen (bingshu.chen@queensu.ca)
}

\note{
This package was build on code developed by Yicong Liu for simple treatment-biomaker interaction model.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bhm}},
\code{\link{coxph}},
\code{\link{lplb.control}}
\code{\link{print.lplb}}
\code{\link{plot.lplb}}
}
\examples{
## fit1 = lplb(Surv(time, status)~z1 + z2 + w, data = dat, B = B, p1 = 2)
## print(fit1)
}
\keyword{ Biomarker interaction, local partial likelihood, bootstrap}
