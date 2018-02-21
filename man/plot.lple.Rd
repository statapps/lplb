\name{plot.lple}
\alias{plot.lple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Plot Function of lple
}
\description{
Draw a series of plots of beta_w vs. w_est for each dependent variable with interactions with the biomarker w. See also:
\code{\link{lple}},
\code{\link{lple.control}}
}
\usage{
\method{plot}{lple}(x, ..., scale = c('original', 'transformed'))
}

\arguments{
  \item{x}{a lple class returned from lple fit.}
  \item{scale}{choose the scale of biomarker variable, 'original' or 'o' for the original biomarker scale. 'transformed' or 't' for transformed scale that mapps biomarker to interval (0, 1). The default is to plot in the original scale.}
  \item{\dots}{other options used in plot().}
}
\details{
plot.lple is called to plot the relationships between beta_w and w_est for each dependent variable with interactions with the biomarker w, from the \code{\link{lple}} fit model.

The number of interaction terms can be set in \code{\link{lple.control}}.

The default method, print.default has its own help page. Use methods("print") to get all the methods for the print generic.
}

\author{Bingshu E. Chen and Siwei Zhang}

\seealso{
\code{\link{lple}},
\code{\link{lple.control}},
\code{\link{print.lple}}
}
\examples{
#
# fit = lplb(y ~ x)
# plot(fit)
#
}

\keyword{plot}
\keyword{lple}
