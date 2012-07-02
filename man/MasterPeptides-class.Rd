\name{MasterPeptides-class}
\Rdversion{1.1}
\docType{class}
\alias{MasterPeptides-class}
\alias{class:MasterPeptides}
\alias{show,MasterPeptides-method}

\title{Class \code{"MasterPeptides"}}

\description{
  A class to store the results of \code{\link{makeMaster}}. This class
  stored the 2 versions (orders) of the \emph{master} final peptide data.
}


\section{Objects from the Class}{
  Objects can be created by calls of the form \code{\link{makeMaster}}.
}

\section{Slots}{
  \describe{
    \item{\code{masters}:}{Object of class \code{"list"} storing the 2
      master \code{data.frame} obejcts. }
    \item{\code{pepfiles}:}{Object of class \code{"character"} with the
      list of final peptide iinput files. }
    \item{\code{orders}:}{Object of class \code{"list"} with the
      \code{numeric} vectors specifying the order of \code{pepfiles}
      used to generate the respective \code{masters}
      \code{data.frames}. }
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "MasterPeptides")}: to print a
      textual representation of the instance. }
  }
}

\references{
  Shliaha P.V., Bond N. J., Lilley K.S. and Gatto L., in prep.
}

\author{
  Laurent Gatto
}


\keyword{classes}
