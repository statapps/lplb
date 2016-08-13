\name{plot.lplb}
\alias{plot.lplb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Plot Function of lplb
}
\description{
Draw a series of plots of beta_w vs. w_est for each dependent variable with interactions with the biomarker w. See also:
\code{\link{lplb}},
\code{\link{lplb.control}}
}
\usage{
\method{plot}{lplb}(x, ...)
}

\arguments{
  \item{x}{a lplb class returned from lplb fit}
  \item{\dots}{other options used in plot()}
}
\details{
plot.lplb is called to plot the relationships between beta_w and w_est for each dependent variable with interactions with the biomarker w, from the \code{\link{lplb}} fit model.

The number of interaction terms can be set in \code{\link{lplb.control}}.

The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}

\author{Bingshu E. Chen and Siwei Zhang}

\seealso{
\code{\link{lplb}},
\code{\link{lplb.control}},
\code{\link{print.lplb}}
}
\examples{
#
# plot.lplb(x)
#
}

\keyword{plot}
\keyword{lplb}
