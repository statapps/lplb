\name{print.lple}
\alias{print.lple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
print a lple object
}
\description{
print are used to provide a short summary of lple outputs.
}
\usage{
\method{print}{lple}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a lple class returned from lple fit}
  \item{\dots}{other options used in print()}
}
\details{
  print.lple is called to print object or summary of object from the lple model \code{\link{lple}}.

  The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}
\author{Siwei Zhand and Bingshu E. Chen}
\seealso{
 The default method for print \code{\link{print.default}},
\code{\link{lple}}
}
\examples{
#
#  print(fit)
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ print }

