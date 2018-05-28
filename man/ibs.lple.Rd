\name{ibs}
\alias{ibs.default}
\alias{ibs.lple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calcuate integrated Brier Score for lple
}
\description{
The function ibs is used to calcualte integrated Brier Score for lple.
}
\usage{
\method{ibs}{lple}(object, newdata=NULL, newy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the results of a lple fit}
  \item{newdata}{If no null, new data and new y will be used for calculation}
  \item{newy}{Used together with new data}
}
\details{
  ibs.lple is called to calculate integrate Brier Score for the lple model \code{\link{lple}}.

  The default method, ibs.default has its own help page. Use methods("ibs") to get all the methods for the ibs generic.
}
\author{Wen Teng and Bingshu E. Chen}
\seealso{
 The default method for ibs \code{\link{ibs.default}},
\code{\link{lple}}
}
\examples{
### To calculate Brier score for the original fitted data
#
#  ibs(fit)
#
### To calculate Brier score for new data with new outcomes
#
#  ibs(fit, newdata, newy)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Integrated Brier Score }

