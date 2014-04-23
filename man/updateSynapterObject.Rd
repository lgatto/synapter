\name{updateSynapterObject}
\alias{updateSynapterObject}
\title{Updates an old synapter object}
\usage{
updateSynapterObject(obj, verbose = TRUE)
}
\arguments{
  \item{obj}{synapter object}

  \item{verbose}{print verbose output}
}
\value{
new synapter object
}
\description{
This function updates an old synapter object. Please see
the details section for more information about the changes
in different versions of synapter objects.
}
\details{
\describe{ \item{synapter 1.7.1}{ Introduce
\code{updateSynapterObject} function. \cr Add
IdentSpectrumFile, IdentSpectrumData, QuantSpectrumFile,
QuantSpectrumData, IdentFragmentFile, IdentFragmentData,
QuantFragmentFile, and QuantFragmentData fields to synapter
class.} }
}
\examples{
synapterTinyData()
synapterTiny
newSynapterObj <- updateSynapterObject(synapterTiny)
newSynapterObj
}
\seealso{
\code{\link{Synapter}}
}

