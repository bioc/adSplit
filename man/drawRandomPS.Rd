\name{drawRandomPS}
\alias{drawRandomPS}
\title{Draw sets of probe-sets}
\description{
  This function draws a given number of probe-sets randomly, such that
  probe-sets referring to the same are either included or excluded as
  a whole.
}
\usage{
drawRandomPS(nps, EID2PSenv, allEIDs)
}
\arguments{
  \item{nps}{number of probe-sets to be drawn.}
  \item{EID2PSenv}{a hash mapping EntrezGene to probe-set identifiers.}
  \item{allEIDs}{vector of all EntrezGene identifiers represented on a chip.}
}
\value{
  A named vector of probe-set identifiers. The names correspond to the
  EntrezGene identifiers.
}
\author{Claudio Lottaz}
\examples{
# draw ten random probe-sets from hu6800
library(hu6800.db)
EID2PSenv <- makeEID2PROBESenv(hu6800ENTREZID)
drawRandomPS(10, EID2PSenv, ls(EID2PSenv))
}
\keyword{datagen}
