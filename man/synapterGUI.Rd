\name{synapterGUI}
\alias{synapterGUI}
\title{Launch the synapter GUI}
\usage{
  synapterGUI(n = 1)
}
\arguments{
  \item{n}{Number of analysis to be run.}
}
\value{
  Invisibly returns \code{null}. Used for its side effects.
}
\description{
  This function starts the synapter graphical user
  interface that allows to (1) load \code{n} sets of files
  to be analyses, (2) parameters used for data filtering
  and retention time modelling and (3) parameters used for
  grid search. The \code{n} file sets and the parameters
  recorded by the GUI are used to call
  \code{\link{synergise}} that performs the analysis and
  generates \code{n} html reports. See the \code{synapter}
  vignette for a description of the interface and
  \code{\link{synergize}} for a description of the
  parameters and data processing.
}
\author{
  Laurent Gatto
}

