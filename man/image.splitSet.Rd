\name{image.splitSet}
\alias{image.splitSet}
\title{Illustrate Split Sets}
\description{
  Draws an image of all splits, one per row, of a \code{splitSet}
  object. Each column corresponds to a patient.
}
\usage{
image.splitSet(x, filter.fdr = 1, main = "", show.graph = FALSE, 
               max.label.length = 50, full.names = TRUE, 
               xlab = NULL, sample.labels = FALSE, 
               col = c("yellow", "red"), invert = FALSE, 
               outfile = NULL, res = 72, pointsize = 7, ...)
}
\arguments{
  \item{x}{the object of class \code{splitSet} to be illustrated.}
  \item{filter.fdr}{worst acceptable false discovery rate for the
    shown set of splits. All splits with q-values beloww this level are
    dropped from the image.}
  \item{main}{a title for the image.}
  \item{show.graph}{whether hierarchical relations between GO terms
    should be illustrated in a graph-like display to the left of the
    image.}
  \item{max.label.length}{Maximal length of the annotations shown to
    the right of the image. Longer annotations are truncated.}
  \item{full.names}{Show full names for annotations instead of their
    identifiers only.}
  \item{xlab}{additional annotation on the x-axis.}
  \item{sample.labels}{whether names of samples are to be shown on the
    x-axis.}
  \item{col}{two strings encoding the colors to be used to illustrate
    to which group a sample is attributed.}
  \item{invert}{whether to draw in white on black background.}
  \item{outfile}{the filename on which to draw the image in postscript
    format. The default is \code{NULL}, meaning to produce the image
    interactively.}
  \item{res}{resolution for bitmap output on postscript.}
  \item{pointsize}{size of font.}
  \item{\dots}{further arguments passed to \code{image}.}
}
\details{
  The set of splits given is illustrated as an image. Each row
  corresponds to an annotation, each column to a patient. In position
  (x,y), the association of patient x to a group with respect to
  annotation y is coded as colors (yellow and red by default). The
  image is ordered by hierarchical clustering such that similar patients
  and similar splits are brought closer together.
}
\value{
  Always returns NULL.
}
\author{Claudio Lottaz}
\seealso{\code{\link{adSplit}}}
\examples{
data(golubKEGGSplits)
image(golubKEGGSplits, filter.fdr=0.5)
}
\keyword{hplot}



