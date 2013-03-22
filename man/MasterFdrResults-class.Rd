\name{MasterFdrResults-class}
\Rdversion{1.1}
\docType{class}

\alias{MasterFdrResults-class}
\alias{class:MasterFdrResults}
\alias{MasterFdrResults}
\alias{bestComb,MasterFdrResults-method}
\alias{bestComb}
\alias{allComb,MasterFdrResults-method}
\alias{allComb}
\alias{masterFdr,MasterFdrResults-method}
\alias{masterFdr}
\alias{fileNames,MasterFdrResults-method}
\alias{plot,MasterFdrResults,missing-method}
\alias{show,MasterFdrResults-method}

\title{Class \code{"MasterFdrResults"}}

\description{
  This class stored the results of the 
}


\section{Objects from the Class}{
Objects are created with the \code{\link{estimateMasterFdr}} function.
}

\section{Slots}{
  \describe{
    \item{\code{all}:}{Object of class \code{"data.frame"} with results
      of all possible combinations. }
    \item{\code{files}:}{Object of class \code{"character"} with file
      names used in combinations. }
    \item{\code{best}:}{Object of class \code{"data.frame"} with results
      of the best combination. }
    \item{\code{masterFdr}:}{ Object of class \code{"numeric"} storing
      the best combination. }
  }
}
\section{Methods}{
  \describe{
    \item{bestComb}{\code{signature(object = "MasterFdrResults")}: ... }
    \item{allComb}{\code{signature(object = "MasterFdrResults")}: ... }
    \item{fileNames}{\code{signature(object = "MasterFdrResults")}: ... }
    \item{masterFdr}{\code{signature(object = "MasterFdrResults")}: ... }
    \item{plot}{\code{signature(x = "MasterFdrResults", y = "missing")}: ... }
    \item{show}{\code{signature(object = "MasterFdrResults")}: ... }
  }
}

\references{
Improving qualitative and quantitative performance for MSE-based label
free proteomics, N.J. Bond, P.V. Shliaha, K.S. Lilley and L. Gatto,
Journal of Proteome Research, 2013, in press.

The Effects of Travelling Wave Ion Mobility Separation on Data
Independent Acquisition in Proteomics Studies, P.V. Shliaha, N.J. Bond,
L. Gatto and K.S. Lilley, Journal of Proteome Research, 2013, in press. 
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\examples{
\dontrun{
library(synapterdata)
loadMaster()
class(master)
master
}
}

\keyword{classes}
