\name{synapterTinyData}
\alias{synapterTinyData}
\alias{synapterTiny}

\title{
  Loads a small test data for the 'synapter' package
}

\description{
  Instead of using \code{data} to load the \code{synapterTiny} data set,
  \code{synapterTinyData} will load it and initialise it for proper
  downstream analysis, during which the \code{04_test_database.fasta}
  file, provided with the package and references inside the object
  needs, to be accessed. However, as the exact location can not be known
  in advance, the reference is updated with the file's correct local
  path. 
  
  This data set has been generated with the \code{\link{Synapter}}
  constructor. Note that the input data file sizes have been reduced by
  depleting many rows (peptides and EMRTs) from the original
  csv files.
  
  In addition, several columns that where not necessary for
  processing were also removed. As such, the data stored in 
  \code{synapterTiny} does not relect the data obtained when following
  the section 'Preparing the input data' in the section, without however
  affecting the processing and final results. 
}

\usage{
synapterTinyData()
}


\value{
  A character vector with the data set name, "synapterTiny". Used
  for its side effect of loading \code{synapterTiny}, an instance of
  class \code{\link{Synapter}}, in \code{.GlovalEnv}.
}

\source{
  Bond, Shliaha, Gatto and Lilley, \emph{in prep.}
}

\author{
Laurent Gatto
}


\examples{
synapterTinyData()
synapterTiny
}

\keyword{datasets}