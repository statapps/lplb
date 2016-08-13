\name{print.lplb}
\alias{print.lplb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a lplb object
}
\description{
print are used to provide a short summary of lplb outputs.
}
\usage{
\method{print}{lplb}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a lplb class returned from lplb fit}
  \item{\dots}{other options used in print()}
}
\details{
  print.lplb is called to print object or summary of object from the lplb model \code{\link{lplb}}.

  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}
\author{Siwei Zhand and Bingshu E. Chen}
\seealso{
 The default method for print \code{\link{print.default}},
\code{\link{lplb}}
}
\examples{
#
#  print(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }

