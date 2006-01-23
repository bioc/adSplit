\name{makeLID2PROBESenv}
\alias{makeLID2PROBESenv}
\title{Generate LID2PROBES environment}
\description{
  Make hash containing probe-sets per LOCUSLINK identifier.
}
\usage{
makeLID2PROBESenv(LIDenv)
}
\arguments{
  \item{LIDenv}{an environment containing one entry per probe-set
    holding all corresponding LOCUSLINK identifiers.}
}
\value{
  An environment containing one entry per LOCUSLINK identifier holding
  all corresponding probe-sets.  
} 
\author{Joern Toedling, Claudio Lottaz} 
\examples{ 
library(hu6800)
makeLID2PROBESenv(hu6800LOCUSID) 
}

\keyword{datagen}
