dbUniquePeptideSet <- function(fastafile, missedCleavages = 0, PLGS = TRUE,
                               verbose = TRUE) {

    missedCleavages <- max(missedCleavages)
    stopifnot(missedCleavages >= 0)

    proteins <- as.character(readAAStringSet(fastafile))

    if (PLGS) {
      ## we modify the standard "cleaver" rules to emulate PLGS' behaviour
      ## according to the discussion on:
      ## https://github.com/sgibb/cleaver/issues/4 and
      ## https://github.com/lgatto/synapter/issues/53
      ##
      ## 1. We cleave after every K or R that is *not* followed by P (default
      ##    trypsin rule), or D, E, K, R (PLGS rule)
      ##    [perfect cleavage, no mis-cleavages allowed].
      ## 2. Subsequently we count the consecutive K or R in each peptide and
      ##    consider these numbers as maximal possible mis-cleavages to get
      ##    every peptide of interest.
      ## 3. We compare the numbers from 2 with the user-defined threshold
      ##    "missedCleavages" and only allow 0 to
      ##    min(number from 2, user-defined threshold).
      ## 4. We cleave after every K or R that is *not* followed by P (default
      ##    trypsin rule) and allow maximal n mis-cleavages. n is determined
      ##    in 3.
      ##
      ## See also:
      ##
      ## Glatter, Timo, et al. Large-scale quantitative assessment of
      ## different in-solution protein digestion protocols reveals superior
      ## cleavage efficiency of tandem Lys-C/trypsin proteolysis over trypsin
      ## digestion.
      ## Journal of proteome research 11.11 (2012): 5145-5156.
      ## http://dx.doi.org/10.1021/pr300273g
      ##
      ## Rodriguez, Jesse, et al. Does trypsin cut before proline?.
      ## Journal of proteome research 7.01 (2007): 300-305.
      ## http://dx.doi.org/10.1021/pr0705035
      ##
      ## Brownridge, Philip, and Robert J. Beynon. The importance of the digest:
      ## proteolysis and absolute quantification in proteomics.
      ## Methods 54.4 (2011): 351-360.
      ## http://dx.doi.org/10.1016/j.ymeth.2011.05.005
      ##

      ## 1) "perfect split" without any mis-cleavage
      peptides <- unlist(cleave(proteins, custom="[K|R](?=[^PDEKR])",
                                missedCleavages=0), use.names=FALSE)

      ## 2) guess maximal number of maximal needed mis-cleavages
      rx <- gregexpr("[KR]+", peptides)
      maxMissedCleavages <- sapply(rx, function(x)max(attr(x, "match.length")))

      ## We reduce maxMissedCleavages by one because we assume that all peptides
      ## that are cleaved in 1) have K, R at their right end (and we don't want to
      ## cleave something like "AAAK" again).
      maxMissedCleavages <- maxMissedCleavages-1

      ## 3) determine maximal allowed mis-cleavages
      maxMissedCleavages <- pmin(maxMissedCleavages, missedCleavages)

      ## only cleave peptides that contain at least one K or R
      toCleave <- maxMissedCleavages > 0

      ## 4) final cleaving
      peptides[toCleave] <- unlist(mapply(function(x, m) {
        cleave(x, custom="[K|R](?=[^P])", missedCleavages=(0:m))
      }, x=peptides[toCleave], m=maxMissedCleavages[toCleave],
      SIMPLIFY=FALSE, USE.NAMES=FALSE), use.names=FALSE)
    } else {
      ## using default cleaver cleavages rules for trypsin, defined:
      ## http://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#Tryps
      peptides <- unlist(cleave(proteins, enzym="trypsin",
                                missedCleavages=0:missedCleavages),
                         use.names=FALSE)
    }

    ## remove non-unique peptides
    upeptides <- peptides[!(duplicated(peptides) |
                            duplicated(peptides, fromLast=TRUE))]

    ## Note that this does however not take pre/post-fix peptides into account
    ## like "DLIELTESLFR" and "HNPEFTMMELYMAYADYHDLIELTESLFR",
    ##      "YYGYTGAFR" and "EGYYGYTGAFR"
    ## where the first ones are NOT unique!
    ## Such cases are handled by filtering duplicates in the peptide data
    if (verbose)
        message(paste("Fasta file: ", length(proteins), " proteins\n",
                      "            ", length(upeptides),
                      " out of ", length(unique(peptides)),
                      " tryptic peptides are proteotypic.\n",
                      "Used rule:  ", ifelse(PLGS, "PLGS", "cleaver"),
                      sep = ""))
    return(upeptides)
}

