\name{createUniquePeptideDbRds}
\alias{createUniquePeptideDbRds}
\title{Create an RDS file for the 'Unique Peptides Database'}
\usage{
createUniquePeptideDbRds(fastaFile, outputFile = paste0(fastaFile, ".rds"),
  missedCleavages = 0, IisL = FALSE, verbose = interactive())
}
\arguments{
  \item{fastaFile}{file path of the input fasta file}

  \item{outputFile}{file path of the target RDS file; must
  have the file extension ".rds"}

  \item{missedCleavages}{Number of maximal allowed missed
  cleavages. Default is 0.}

  \item{IisL}{If \code{TRUE} Isoleucin and Leucin are
  treated as identical. In this case sequences like "ABCI",
  "ABCL" are removed because they are not unqiue. If
  \code{FALSE} (default) "ABCI" and "ABCL" are reported as
  unique.}

  \item{verbose}{If \code{TRUE} a verbose output is
  provied.}
}
\description{
This function creates an RDS file to store the
\dQuote{Unique Peptides Database} for future runs of
\code{\link{synergise}} or \code{\link{Synapter}}.
}
\examples{
\dontrun{
createUniquePeptideDbRds("uniprot.fasta", "uniprot.fasta.rds")
}
}
\author{
Sebastian Gibb <mail@sebastiangibb.de>
}
\seealso{
\code{\link{Synapter}} for details about the cleavage
procedure.
}

