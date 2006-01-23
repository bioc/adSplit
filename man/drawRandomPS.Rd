\name{drawRandomPS}
\alias{drawRandomPS}
\title{Draw sets of probe-sets}
\description{
  This function draws a given number of probe-sets randomly, such that
  probe-sets referring to the same are either included or excluded as
  a whole.
}
\usage{
drawRandomPS(nps, LID2PSenv, allLIDs)
}
\arguments{
  \item{nps}{number of probe-sets to be drawn.}
  \item{LID2PSenv}{a hash mapping LOCUSLINK to probe-set identifiers.}
  \item{allLIDs}{vector of all LOCUSLINK identifiers represented on a chip.}
}
\value{
  A named vector of probe-set identifiers. The names correspond to the
  LOCUSLINK identifiers.
}
\author{Claudio Lottaz}
\examples{
# draw ten random probe-sets from hu6800
library(hu6800)
LID2PSenv <- makeLID2PROBESenv(hu6800LOCUSID)
drawRandomPS(10, LID2PSenv, ls(LID2PSenv))
}
\keyword{datagen}
