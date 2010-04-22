\name{diana2means}
\alias{diana2means}
\title{2-Means with Hierarchical Initialization}
\description{
  Split a set of data points into two coherent groups using the
  k-means algorithm. Instead of random initialization, divisive
  hierarchical clustering is used to determine initial groups and the
  corresponding centroids.
}
\usage{
diana2means(mydata, mingroupsize = 5, 
            ngenes = 50, ignore.genes = 5, 
            return.cut = FALSE)
}
\arguments{
  \item{mydata}{either an expression set as defined by the package
    \code{Biobase} or a matrix of expression levels (rows=genes,
    columns=samples).}
  \item{mingroupsize}{report only splits where both groups are larger than this size.}
  \item{ngenes}{number of genes used to compute cluster quality DLD-score.}
  \item{ignore.genes}{number of best scoring genes to be ignored when
     computing DLD-scores.}
  \item{return.cut}{logical, whether to return the attributions of samples
     to groups.}  
} 
\details{

  This function uses divisive hierarchical clustering (diana) to
  generate a first split of the data. Thereby, each column of the data
  matrix is considered to represent a data element. From the thus
  generated temptative groups, centroids are deduced and used to
  initialize the k-means clustering algorithm.

  For the split optimized by k-means the DLD-score is determined using
  the \code{ngenes} and \code{ignore.genes} arguments. 
} 
\value{
  If the logical \code{return.cut} is set to \code{FALSE} (the
  default), a single number is representing the DLD-score for the
  generated split is returned. Otherwise an object of class
  \code{split} containing the following elements is returned:
  \item{cut}{one number out of 0 and 1 per column in the original
    data, specifying the split attribution.}
  \item{score}{the DLD-score achieved by the split.}
}
\author{Joern Toedling, Claudio Lottaz}
\seealso{\code{\link[cluster]{diana}}}
\examples{
# get golub data
library(vsn)
library(golubEsets)
data(Golub_Merge)

# use 10% most variable genes
e <- exprs(Golub_Merge)
vars <- apply(e, 1, var)
e <- e[vars > quantile(vars,0.9),]

# use diana2means to get splits and scores
diana2means(e)
diana2means(e, return.cut=TRUE)
}
\keyword{datagen}

