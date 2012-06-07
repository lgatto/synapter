\name{synapter-package}
\alias{synapter-package}
\alias{synapter}
\docType{package}

\title{
  Combines label-free data for optimal identification and quantitation 
}

\description{
  The synapter package provides functionality to reanalyse proteomics
  data acquired on a Synapt G2 mass spectrometer. One or several 
  files possibly run with additional ion mobility separation can be
  combined to other files to optimise identification and quantitation.  
}


\details{
  Three pipelines of variying flexibility are proposed the preform data
  analysis: (1) a graphical user interface is provided with the
  \code{synatperGUI} function, (2) the \code{synergise} function is a
  single entry function for a complete analysis and (3) low level,
  step-by-step data processing can be achieved as described in 
  \code{?Synapter}. 
  
  A high-level overview of the package and how to operate it can be
  found in the vignette, accessinble with
  \code{synapterGuide()}. Detailed information about the data processing 
  can be found in the respective function and class manual pages
  appropriately referenced in the vignette. 
}

\author{
Laurent Gatto, Pavel V. Shliaha and Nick J. Bond

Maintainer: Laurent Gatto <lg390@cam.ac.uk>
}

\keyword{ package }

