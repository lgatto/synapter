#' create fragment library out of multiple fragment files
#'
#' @param master master data.frame
#' @param files path to final_fragment files
#' @param removeNeutralLoss remove rows with neutral loss != "none"?
#' @param verbose verbose output?
#'
#' @return a fragment library data.frame (similar to the input files but
#' combined)
#' @noRd
.createFragmentLibrary <- function(master, files, removeNeutralLoss=FALSE,
                                   verbose=interactive()) {
  fragments <- lapply(files, function(f) {
    r <- .readFragments(file=f, removeNeutralLoss=removeNeutralLoss,
                        verbose=verbose)
    ## keep only peptides we have seen in master
    r <- r[r$peptide.seq %in% master$peptide.seq, ]
    r
  })

  ## find individual runs
  runs <- rep.int(1:length(files), sapply(fragments, nrow))

  fragments <- do.call(rbind, fragments)

  ## store individual runs in df
  fragments$run <- runs

  ## recalculate intensities
  if (verbose) {
    message("Recalculate fragments' intensity values.")
  }
  fragments <- .recalculateFragmentIntensities(fragments)

  ## match precursor.leID
  if (verbose) {
    message("Regenerate precursor.leID values.")
  }
  i <- match(fragments$peptide.seq, master$peptide.seq)
  fragments$precursor.leID <- master$precursor.leID[i]

  rownames(fragments) <- NULL

  fragments
}

#' recalculate fragment intensities for a library
#'
#' @param fragments data.frame read by .createFragmentLibrary
#'
#' @return a data.frame with recalculated intensity values
#' @noRd
.recalculateFragmentIntensities <- function(fragments) {
  ## sum fragment intensity for equal fragments of each peptide
  product.inten <- ave(fragments$product.inten,
                       fragments$peptide.seq, fragments$fragment.str,
                       FUN=sum)
  ## sum precursor intensity for equal fragments of each peptide
  precursor.inten <- ave(fragments$precursor.inten,
                         fragments$peptide.seq, fragments$fragment.str,
                         FUN=sum)
  ## calculate new intensity
  fragments$product.inten <- product.inten / precursor.inten

  ## calculate new precursor intensity for each peptide
  split(fragments$precursor.inten, fragments$peptide.seq) <-
    lapply(split(fragments[c("precursor.inten", "run")], fragments$peptide.seq), function(x) {
      round(mean(x$precursor.inten[!duplicated(x$run)]))
  })

  ## take only the first occurence of each fragment
  fragments <-
    fragments[!duplicated(paste0(fragments$peptide.seq, fragments$fragment.str)), ]

  ## normalize product intensity by new precursor intensity
  ## and round to nearest integer
  fragments$product.inten <- round(fragments$precursor.inten * fragments$product.inten)
  fragments$product.rank <- ave(fragments$product.rank,
                                fragments$peptide.seq,
                                FUN=seq_along)
  rownames(fragments) <- NULL
  fragments
}
