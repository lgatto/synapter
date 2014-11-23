\name{plotCrossMatching}
\alias{plotCrossMatching}
\alias{plotCrossMatching,Synapter-method}
\docType{methods}

\title{Plot cross matching results.}

\description{
  This method plot the results of the cross matching procedure
  (\code{\link{crossMatching}}). A single plot contains two pnaels. The upper
  panel shows the identification fragments and the lower one the MS2 spectrum of
  the quantitation run. Common peaks are drawn in a slightly darker colour and
  with a point on the top.
}

\arguments{
  \item{object}{ Object of class \code{"\link{Synapter}"} .}
  \item{key}{ \code{character}, value to look for. }
  \item{column}{ \code{character}, name of the column in which
    \code{plotCrossMatching} looks for \code{key}.}
  \item{\dots}{Further arguments passed to internal functions.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "Synapter", key = "character",
        column = "character", verbose = "logical", \dots)}}{
      Plots identification fragments against quantitation spectra.
      The \code{\dots} arguments are passed to the internal functions.
      Currently \code{legend.cex}, \code{fragment.cex}, and
      most of the \code{\link[graphics]{plot.default}} arguments
      (like \code{xlim}, \code{ylim}, \code{main}, \code{xlab}, \code{ylab},
      \dots) are supported.
      \code{legend.cex} and \code{fragment.cex} control the size of the legend
      and fragments labels (default: \code{0.5}).
      Please see \code{\link[graphics]{par}} for details about \code{cex}.
      If \code{verbose = TRUE} a progress bar is shown.
    }
  }
}

\seealso{
  \code{\link{crossMatching}}
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\keyword{methods}

