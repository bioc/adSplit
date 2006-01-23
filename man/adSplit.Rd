\name{adSplit}
\alias{adSplit}
\title{Annotation-Driven Splits}
\description{
  This function searches for annotation-driven splits of patients in
  microarray data. A split is a partitionning of patients into two
  groups. In order to do so it referes to GO terms and KEGG
  pathways. In addition, a significance measure can be computed by
  simulating a random distribution of scores. DLD-scores are used 
  to judge the quality of a split.
}
\usage{
adSplit(mydata, annotation.ids, chip.name, 
        min.probes = 20, max.probes = NULL, 
        B = NULL, min.group.size = 5, ngenes = 50, 
        ignore.genes = 5)
}
\arguments{
  \item{mydata}{either an expression set as defined by the package
    \code{Biobase} or a matrix of expression levels (rows=genes,
    columns=samples).}
  \item{annotation.ids}{a vector of GO or KEGG identifiers in the form
    "GO:..." or "KEGG:..." respectively. The prefix "KEGG:" is removed
    from the KEGG-identifiers before accessing the chip's "...PATH2PROBES"
    hash.}
  \item{chip.name}{the name of the chip by which the expression set is
    measured. \code{adSplit} attempts to load a library of the same name
    and expects to find a hash called "<chip-name>GO2ALLPROBES" and one
    called "<chip-name>PATH2PROBES" there.}
  \item{min.probes}{annotation identifiers with fewer than this
    associated genes are skipped.}
  \item{max.probes}{annotation identifiers with more than this
    associated genes are skipped. The default is ten percent of the
    genes on the chip.}
  \item{B}{the number of random gene set samplings to be performed to 
    compute empirical p-values.}
  \item{min.group.size}{filter criteria to avoid splits suggesting
    tiny groups. Splits where one of the two suggested groups are smaller
    than this number are removed from the split set.}
  \item{ngenes}{number of genes used to compute DLD scores.}
  \item{ignore.genes}{number of best scoring genes to be ignored when
    computing DLD scores.} 
}
\details{
  This function applies the same splitting procedure to all annotation
  identifiers provided. Firstly, the associated genes for one identifier
  are determined and extracted from the expression data. Then the
  \code{diana2means} function is applied to the restricted data and the
  different splits generated are collected into a single \code{splitSet}
  object.

  As annotation identifiers vectors of identifiers of the
  \code{KEGG:nnnnn} and \code{GO:nnnnnn} are valid. In addition, the
  keywords "KEGG", "GO" and "all" are allowed, representing all terms in
  the corresponding ontology. 

  If \code{B} is set to a integer number this number of samplings are
  used to generate a null-distribution of DLD-scores. This
  distribution is used to compute empirical p-values for each
  split. If more than one valid split is found, multiple testing is
  corrected for by applying Benjamini-Hochbergs correction from the
  multtest package.
} 
\value{
  Returns an object of class \code{splitSet} with the following list
  elements: 
  \item{cuts}{a matrix of split attributions. One row per annotation
    identifier (GO term or KEGG pathway for which a split has been
    generated. One column per object in the dataset.}
  \item{score}{one score per generated split.}
  \item{pvalue}{one empirical p-value per generated split, or \code{NULL}}
  \item{qvalue}{one q-value computed according Benjamini-Hochberg's
    correction for multiple testing per generated split, or \code{NULL}}
}
\author{Claudio Lottaz, Joern Toedling}
\seealso{\code{\link{diana2means}}, \code{\link{randomDiana2means}}, 
         \code{\link{image.splitSet}}}
\examples{ 
# prepare data
library(vsn) 
library(golubEsets) 
data(golubMerge) 
golubNorm <- vsn(golubMerge, describe.preprocessing=FALSE, verbose=FALSE)

# generate annotation-driven splits for apoptosis and signal transduction
x <- adSplit(golubNorm, "GO:0006915", "hu6800")
x <- adSplit(golubNorm, c("GO:0007165","GO:0006915"), "hu6800", max.probes=7000)

# generate a split for glutamate metabolism including 
# an empirical p-value
x <- adSplit(golubNorm, "KEGG:00251", "hu6800", B=100)

# generate splits for all KEGG pathways.
x <- adSplit(golubNorm, "KEGG", "hu6800")
image(x)
}
\keyword{datagen}

