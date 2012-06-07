readFasta <- function(file, checkComments = TRUE)  {
  ## Adapted from Biostrings::readFASTA (deprecated anyway)
  if (!file.exists(file))
    stop(file, " file not found.")
  s1 <- scan(file = file, what = "", sep = "\n", quote = "", 
             allowEscapes = FALSE, quiet = TRUE)
  if (checkComments) {
    comments <- grep("^;", s1)
    if (length(comments) > 0) 
      s1 <- s1[-comments]
  }
  descriptions <- which(substr(s1, 1L, 1L) == ">")
  numF <- length(descriptions)
  if (numF == 0) 
    stop("no FASTA sequences found")
  dp <- descriptions + 1L
  dm <- descriptions - 1L
  end <- c(dm[-1], length(s1))  
  lapply(seq_len(numF), function(i) {
        desc <- s1[descriptions[i]]
        desc <- substr(desc, 2L, nchar(desc)) ## stripping leading '^>'
        if (end[i] >= dp[i]) {
            seq <- paste(s1[dp[i]:end[i]], collapse = "")
        }
        else {
            warning("record \"", desc, "\" contains no sequence")
            seq <- ""
        }
        list(desc = desc, seq = seq)
    })
}


digest <- function (sequence,
                    enzyme = "trypsin", ## curently ignored
                    missed = 0) {       ## currently ignored  
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  enzyme <- "trypsin"
  if (seq_vector[end_position] == "K" | seq_vector[end_position] == "R") {
    seq_vector[end_position] <- "!"
    seq_string <- paste(seq_vector, collapse = "")
  }
  else seq_string <- sequence
  seq_string <- gsub("KP", "!P", seq_string)
  seq_string <- gsub("RP", "!P", seq_string)
  seq_vector <- strsplit(seq_string, split = "")[[1]]
  stop <- grep("K|R", seq_vector)
  start <- stop + 1
  if (length(stop) == 0) {
    ## warning("sequence does not contain cleavage sites")
    peptides <- character()
  } else {  
    start <- c(1, start)
    stop <- c(stop, end_position)
    peptides <- substring(sequence, start, stop)
  }
  return(peptides)
}

dbUniquePeptideSet <- function(fastafile, verbose = TRUE) {
  db <- readFasta(fastafile)
  seq <- sapply(db, function(x) x$seq)
  peptides <- unlist(sapply(seq, digest))
  names(peptides) <- NULL
  tab <- table(peptides)
  upeptides <- names(tab[tab == 1])
  ## Note that this does however not take pre/post-fix peptides into account
  ## like "DLIELTESLFR" and "HNPEFTMMELYMAYADYHDLIELTESLFR",
  ##      "YYGYTGAFR" and "EGYYGYTGAFR"  
  ## where the first ones are NOT unique!
  ## Such cases are handled by filtering duplicates in the peptide data
  if (verbose)
    message(paste("Fasta file: ", length(seq), " proteins\n", 
                  "            ", length(upeptides),
                  " out of ", length(peptides),
                  " tryptic peptides are unique.",
                  sep = ""))
  return(upeptides)
}
  
