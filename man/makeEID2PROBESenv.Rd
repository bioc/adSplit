\name{makeEID2PROBESenv}
\alias{makeEID2PROBESenv}
\title{Generate EID2PROBES environment}
\description{
  Make hash containing probe-sets per EntrezGene identifier.
}
\usage{
makeEID2PROBESenv(EIDenv)
}
\arguments{
  \item{EIDenv}{an environment containing one entry per probe-set
    holding all corresponding EntrezGene identifiers.}
}
\value{
  An environment containing one entry per EntrezGene identifier holding
  all corresponding probe-sets.  
} 
\author{Joern Toedling, Claudio Lottaz} 
\examples{ 
library(hu6800)
makeEID2PROBESenv(hu6800ENTREZID) 
}

\keyword{datagen}
