dbUniquePeptideSet <- function(fastafile, verbose = TRUE) {
  proteins <- as.character(readAAStringSet(fastafile))
  peptides <- unlist(cleave(proteins, enzym="trypsin", missedCleavages=0),
                     use.names=FALSE)
  upeptides <- peptides[!(duplicated(peptides) | duplicated(peptides, fromLast=TRUE))]
  ## Note that this does however not take pre/post-fix peptides into account
  ## like "DLIELTESLFR" and "HNPEFTMMELYMAYADYHDLIELTESLFR",
  ##      "YYGYTGAFR" and "EGYYGYTGAFR"  
  ## where the first ones are NOT unique!
  ## Such cases are handled by filtering duplicates in the peptide data
  if (verbose)
    message(paste("Fasta file: ", length(proteins), " proteins\n", 
                  "            ", length(upeptides),
                  " out of ", length(unique(peptides)),
                  " tryptic peptides are proteotypic.",
                  sep = ""))
  return(upeptides)
}
   
