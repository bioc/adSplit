\name{hist.splitSet}
\alias{hist.splitSet}
\title{Overview Histogram for splitSets}
\description{
  Draws a histogram of empirical p-values and shows the corresponding
  q-values corrected for multiple testing.
}
\usage{
hist.splitSet(x, main = "Distribution of p-Values", 
              xlab = "p-values", col = "grey", xlim = c(0, 1), ...)
}
\arguments{
  \item{x}{object of type \code{splitSet}. Should hold a considerable
    number of splits.}
  \item{main}{main title of the histogram.}
  \item{xlab}{legend for the x-axis.}
  \item{col}{color for the histogram bars.}
  \item{xlim}{limits for the x-axis (p-values).}
  \item{\dots}{further parameters passed on to the default \code{hist}
    function.}
}
\details{
  This function draws a regular histogram of empirical p-values
  observed in the splitSet at hand. The corresponding q-values,
  corrected by the method suggested by Benjamini-Hochberg, are plottet
  into the same graph. The scale for the q-values is shown at the left
  hand side of the plot.
}
\author{Claudio Lottaz}
\seealso{\code{\link{adSplit}}}
\examples{
data(golubKEGGSplits)
hist(golubKEGGSplits, col="red")
}
\keyword{hplot}
