\name{randomDiana2means}
\alias{randomDiana2means}
\title{Generate null-distributions of DLD-scores}
\description{
  Draws a number of random sets of probe-sets consisting of the needed
  size and applies \code{diana2means} to compute DLD scores.
}
\usage{
randomDiana2means(nprobes, data, chip, ndraws = 10000, 
                  ngenes = 50, ignore.genes = 5)
}
\arguments{
  \item{nprobes}{the size of gene sets.}
  \item{data}{a matrix of expression data, rows corerspond to genes, columns to samples.}
  \item{chip}{the name of the used chip.}
  \item{ndraws}{the number of DLD scores computed.}
  \item{ngenes}{the number of genes used to compute DLD scores (passed
    to \code{diana2means}).}
  \item{ignore.genes}{the number of best scoring genes to be ignored
    when computing DLD scores (passed to \code{diana2means})} 
} 
\details{
  This function uses \code{drawRandomPS} to draw \code{ndraws} gene
  sets. On these it applies \code{diana2means} to determine a
  null-distribution of DLD-scores.
}
\value{
  A vector of DLD-scores.  
}
\author{Joern Toedling, Claudio Lottaz}
\seealso{\code{\link{drawRandomPS}}, \code{\link{diana2means}}}
\examples{
# prepare data
library(vsn) 
library(golubEsets) 
data(golubMerge) 
golubNorm <- vsn(golubMerge, describe.preprocessing=FALSE, verbose=FALSE)

# generate DLD scores
scores <- randomDiana2means(20, exprs(golubNorm), "hu6800", ndraws = 500) 
}
\keyword{datagen}
