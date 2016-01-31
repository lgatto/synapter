##' IMPORTANT: update the initialization value of the ClassVersion field
##' everytime you change this file! (regardless if you just change the code in
##' a method or the complete API).
.Synapter <-
    setRefClass("Synapter",
                fields = list(
                    ## introspection
                    Version = "character",
                    ClassVersion = "character",
                    SynapterLog = "character",
                    Master = "logical",
                    used2 = "logical",
                    ## input data
                    QuantPeptideFile = "character",
                    QuantPeptideData = "data.frame",
                    IdentPeptideFile = "character",
                    IdentPeptideData = "data.frame",
                    QuantPep3DFile = "character",
                    QuantPep3DData = "data.frame",
                    ## fasta db
                    DbFastaFile = "character",
                    ## peptide score data - retained for ploting
                    .QuantPeptideScores = "data.frame",
                    .IdentPeptideScores = "data.frame",
                    ## filtering thresholds
                    PepScoreFdr = "numeric",
                    ProtFpr = "numeric",
                    IdentPpmError = "numeric",
                    QuantPpmError = "numeric",
                    ## merging
                    MergedFeatures = "data.frame",
                    ## rt model
                    LowessSpan = "numeric",
                    RtModel = "list",
                    ## grid search
                    Grid = "list", ## was a "matrix",
                    GridDetails = "list",
                    ## matching EMRTs
                    PpmError = "numeric",
                    RtNsd = "numeric",
                    ImDiff = "numeric",
                    MatchedEMRTs = "data.frame",
                    ## FragmentMatching
                    QuantSpectrumFile = "character",
                    QuantSpectrumData = "MSnExp",
                    IdentFragmentFile = "character",
                    IdentFragmentData = "MSnExp",
                    FragmentMatching = "data.frame",
                    FragmentMatchingPpmTolerance = "numeric"),
                methods = list(
                    initialize = function() {
                        ## IMPORTANT: always increase the ClassVersion field if
                        ## you change anything in this file!
                        .self$ClassVersion <- "2.0.0"

                        .self$Version <- as.character(packageVersion("synapter"))
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Instance created on ", date(), sep=""))
                        .self$used2 <- FALSE
                    },
                    loadFiles = function() {
                        'Read input data files.'
                        .self$IdentPeptideFile <-
                            tk_choose.files("",
                                            caption = "Select identification final peptide file",
                                            multi = FALSE)
                        .self$QuantPeptideFile <-
                            tk_choose.files("",
                                            caption = "Select quantitation final peptide file",
                                            multi = FALSE)
                        .self$QuantPep3DFile <-
                            tk_choose.files("",
                                            caption = "Select quantitation Pep3D file",
                                            multi = FALSE)
                        .self$DbFastaFile <-
                            tk_choose.files("",
                                            caption = "Select fasta file",
                                            multi = FALSE)
                        .self$IdentFragmentFile <-
                            tk_choose.files("",
                                            caption = "Select identification fragments file",
                                            multi = FALSE)
                        .self$QuantSpectrumFile <-
                            tk_choose.files("",
                                            caption = "Select quantitation spectra file",
                                            multi = FALSE)
                    },
                    loadMasterData = function() {
                        message("Reading quantitation final peptide file...")
                        .self$QuantPeptideData <- readFinalPeptides(.self$QuantPeptideFile)
                        .self$QuantPeptideData$errorppm <-
                            error.ppm(obs = .self$QuantPeptideData$precursor.mhp,
                                      theo = .self$QuantPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation peptide data [",
                                                     paste(dim(.self$QuantPeptideData), collapse = ","),
                                                     "]", sep=""))

                        if (length(.self$QuantSpectrumFile)) {
                          loadQuantitationSpectra(.self$QuantSpectrumFile)
                        }
                        if (!.self$Master)
                            stop("Identification final peptide is not a master file")
                        message("Reading master identification peptide file...")
                        ext <- file_ext(.self$IdentPeptideFile)
                        if (tolower(ext) == "csv") {
                            .self$IdentPeptideData <- readFinalPeptides(.self$IdentPeptideFile)
                            .self$SynapterLog <- c(.self$SynapterLog,
                                                   paste("Read master identification peptide (csv) data [",
                                                         paste(dim(.self$IdentPeptideData), collapse = ","),
                                                         "]", sep=""))
                            if (length(.self$IdentFragmentFile)) {
                              loadIdentificationFragments(.self$IdentFragmentFile)
                            }
                        } else if (tolower(ext) == "rds") {
                            masterpeps <- readRDS(.self$IdentPeptideFile)
                            ## testing if masterpeps@masters[[1]] or [[2]] to be used
                            .self$IdentPeptideData <- masterpeps@masters[[1]]
                            .self$PepScoreFdr <- masterpeps@fdr
                            .QuantPeptideData <- .self$QuantPeptideData
                            .self$addQuantIdStats(log = FALSE) ## Quant only
                            .self$filterQuantPepScore(fdrMethod = masterpeps@method, log = FALSE)
                            ## Prot FPR filtering - depending on loadIdentOnly and makeMaster
                            ## PPM error filtering - depending on loadIdentOnly and makeMaster
                            .MergedFeatures <- merge(.self$IdentPeptideData,
                                                     .self$QuantPeptideData,
                                                     by.x = "peptide.seq",
                                                     by.y = "peptide.seq",
                                                     suffixes = c(".ident", ".quant"))
                            .deltaRt <-
                                .MergedFeatures$precursor.retT.ident -
                                    .MergedFeatures$precursor.retT.quant
                            if (all(.deltaRt == 0)) {
                                .self$IdentPeptideData <- masterpeps@masters[[2]]
                                .self$used2 <- TRUE
                            }
                            .self$SynapterLog <- c(.self$SynapterLog,
                                                   paste0("Read master (", ifelse(.self$used2, 2, 1) ,
                                                          ") identification peptide (", ext, ") data [",
                                                          paste(dim(.self$IdentPeptideData), collapse = ","), "]"))
                            .self$PepScoreFdr <- numeric() ## re-initialise
                            .self$QuantPeptideData <- .QuantPeptideData

                            if (length(masterpeps@fragmentlibrary)) {
                              .self$IdentFragmentFile <- .self$IdentPeptideFile
                              .self$IdentFragmentData <- masterpeps@fragmentlibrary
                              .self$SynapterLog <-
                                c(.self$SynapterLog,
                                  paste0("Read master identification fragment data [",
                                         length(.self$IdentFragmentData), "]"))

                            } else if (!length(masterpeps@fragmentlibrary) &&
                                       length(.self$IdentFragmentFile)) {
                              loadIdentificationFragments(.self$IdentFragmentFile)
                            }
                        } else {
                            stop("Master peptide file extention not recognised. Must be 'csv' or 'rds'.")
                        }
                        message("Reading quantitation Pep3D file...")
                        .self$QuantPep3DData <- readPep3D(.self$QuantPep3DFile)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation Pep3D data [",
                                                     paste(dim(.self$QuantPep3DData), collapse = ","),
                                                     "]", sep=""))

                        ## test for correct Pep3D file (closes #42)
                        if (!isCorrespondingPep3DataFile(.self$QuantPeptideData, .self$QuantPep3DData)) {
                          stop("The Pep3D file ", sQuote(.self$QuantPep3DFile),
                               " does not correspond to the given Quantitation Final Peptide file ",
                               sQuote(.self$QuantPeptideFile), "!")
                        }

                        ## getting peptide scores
                        .self$.QuantPeptideScores <-
                            filterPeptideMatchType(.self$QuantPeptideData[, c("peptide.seq",
                                                                              "peptide.score",
                                                                              "peptide.matchType",
                                                                              "protein.dataBaseType")])
                        .self$.QuantPeptideScores <-
                            .self$.QuantPeptideScores[!duplicated(.self$.QuantPeptideScores$peptide.seq), ]
                        .self$.IdentPeptideScores <- ## aleady filtered for matchType
                            .self$IdentPeptideData[, c("peptide.seq",
                                                       "peptide.score",
                                                       "peptide.matchType",
                                                       "protein.dataBaseType")]
                        .self$.IdentPeptideScores <-
                            .self$.IdentPeptideScores[!duplicated(.self$.IdentPeptideScores$peptide.seq), ]
                        ## Housekeeping - Quant only
                        .self$QuantPeptideData$protein.falsePositiveRate <-
                            .self$QuantPeptideData$protein.falsePositiveRate / 100
                        message("Filtering...")
                        ## affects #42
                        .self$filterMismatchingQuantIntensities()
                        .self$filterQuantMatchType() ## Quant only
                        .self$filterQuantFunction()
                        ## must be before duplicated ids are removed and after
                        ## Function filtering
                        .self$QuantPep3DData$isotopicDistr <- .catIsotopeDistr(.self$QuantPep3DData)
                        .self$filterDuplicatedQuantSpectrumIds()
                        message("Computing quantitation identification statistics...")
                        .self$addQuantIdStats() ## Quant only
                    },
                    loadData = function() {
                        if (.self$Master)
                            stop("Identification final peptide is a master file")
                        message("Reading identification final peptide file...")
                        .self$IdentPeptideData <- readFinalPeptides(.self$IdentPeptideFile)
                        .self$IdentPeptideData$errorppm <-
                            error.ppm(obs = .self$IdentPeptideData$precursor.mhp,
                                      theo = .self$IdentPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read identification peptide data [",
                                                     paste(dim(.self$IdentPeptideData), collapse = ","),
                                                     "]", sep=""))
                        message("Reading quantitation final peptide file...")
                        .self$QuantPeptideData <- readFinalPeptides(.self$QuantPeptideFile)
                        .self$QuantPeptideData$errorppm <-
                            error.ppm(obs = .self$QuantPeptideData$precursor.mhp,
                                      theo = .self$QuantPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation peptide data [",
                                                     paste(dim(.self$QuantPeptideData), collapse = ","),
                                                     "]", sep=""))
                        message("Reading quantitation Pep3D file...")
                        .self$QuantPep3DData <- readPep3D(.self$QuantPep3DFile)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation Pep3D data [",
                                                     paste(dim(.self$QuantPep3DData), collapse = ","),
                                                     "]", sep=""))

                        ## test for correct Pep3D file (closes #42)
                        if (!isCorrespondingPep3DataFile(.self$QuantPeptideData, .self$QuantPep3DData)) {
                          stop("The Pep3D file ", sQuote(.self$QuantPep3DFile),
                               " does not correspond to the given Quantitation Final Peptide file ",
                               sQuote(.self$QuantPeptideFile), "!")
                        }

                        ## getting peptide scores
                        .self$.IdentPeptideScores <-
                            filterPeptideMatchType(.self$IdentPeptideData[,c("peptide.seq",
                                                                             "peptide.score",
                                                                             "peptide.matchType",
                                                                             "protein.dataBaseType")])
                        .self$.IdentPeptideScores <-
                            .self$.IdentPeptideScores[!duplicated(.self$.IdentPeptideScores$peptide.seq), ]
                        .self$.QuantPeptideScores <-
                            filterPeptideMatchType(.self$QuantPeptideData[,c("peptide.seq",
                                                                             "peptide.score",
                                                                             "peptide.matchType",
                                                                             "protein.dataBaseType")])
                        .self$.QuantPeptideScores <-
                            .self$.QuantPeptideScores[!duplicated(.self$.QuantPeptideScores$peptide.seq), ]
                        ## Housekeeping
                        .self$QuantPeptideData$protein.falsePositiveRate <-
                            .self$QuantPeptideData$protein.falsePositiveRate / 100
                        .self$IdentPeptideData$protein.falsePositiveRate <-
                            .self$IdentPeptideData$protein.falsePositiveRate / 100
                        message("Filtering...")
                        ## affects #42
                        .self$filterMismatchingQuantIntensities()
                        .self$filterMatchType()
                        .self$filterQuantFunction()
                        ## must be before duplicated ids are removed and after
                        ## Function filtering
                        .self$QuantPep3DData$isotopicDistr <- .catIsotopeDistr(.self$QuantPep3DData)
                        .self$filterDuplicatedQuantSpectrumIds()
                        message("Computing identification statistics...")
                        .self$addIdStats()

                        if (length(.self$IdentFragmentFile)) {
                          loadIdentificationFragments(.self$IdentFragmentFile)
                        }
                        if (length(.self$QuantSpectrumFile)) {
                          loadQuantitationSpectra(.self$QuantSpectrumFile)
                        }
                    },
                    loadIdentificationFragments = function(filename,
                                                           removeNeutralLoss=TRUE,
                                                           removePrecursor=TRUE,
                                                           tolerance=25e-6,
                                                           verbose=TRUE) {
                      .self$IdentFragmentFile <- filename

                      if (length(filename)) {
                        .self$IdentFragmentData <-
                          .finalFragment2spectra(df=.self$IdentPeptideData,
                                                 file=.self$IdentFragmentFile,
                                                 storeAll=FALSE,
                                                 removeNeutralLoss=removeNeutralLoss,
                                                 removePrecursor=removePrecursor,
                                                 tolerance=tolerance,
                                                 verbose=verbose)
                        .self$SynapterLog <-
                          c(.self$SynapterLog,
                            paste0("Read identification fragment data [",
                                   length(.self$IdentFragmentData), "]"))
                      }
                    },
                    loadQuantitationSpectra = function(filename,
                                                       removePrecursor=TRUE,
                                                       tolerance=25e-6,
                                                       verbose=TRUE) {
                      .self$QuantSpectrumFile <- filename

                      if (length(filename)) {
                        .self$QuantSpectrumData <-
                          .spectrumXml2spectra(df=.self$QuantPeptideData,
                                               file=.self$QuantSpectrumFile,
                                               storeAll=TRUE,
                                               removePrecursor=removePrecursor,
                                               tolerance=tolerance,
                                               verbose=verbose)
                        .self$SynapterLog <-
                          c(.self$SynapterLog,
                            paste0("Read quantitation spectra [",
                                   length(.self$QuantSpectrumData), "]"))
                      }
                    },
                    getMaster = function() {
                        ' Gets Master field.'
                        return(.self$Master)
                    },
                    setMaster = function(master) {
                        ' Sets Master field.'
                        .self$Master <- master
                        if (master)
                            .self$SynapterLog <- c(.self$SynapterLog,
                                                   "Identification final peptide is a _master_ file.")
                    },
                    addIdStats = function() {
                        'Add p-values and q-values to final peptide data sets.'
                        .self$addQuantIdStats()
                        .self$addIdentIdStats()
                    },
                    addQuantIdStats = function(log = TRUE) {
                        'Add p-values and q-values to final quantitation peptide data sets.'
                        quantStats <- try(getIdStats(.self$QuantPeptideData), silent=TRUE)
                        if (inherits(quantStats, "try-error")) {
                            warning("No Regular and Random protein database types in quantitation peptide data.")
                        } else {
                            .self$QuantPeptideData <- cbind(.self$QuantPeptideData, quantStats)
                            if (log)
                                .self$SynapterLog <- c(.self$SynapterLog,
                                                       "Added identification statistics to quantitation data")
                            .self$filterQuantRandomEntries()
                        }
                    },
                    addIdentIdStats = function() {
                        'Add p-values and q-values to identification final peptide data sets.'
                        identStats <- try(getIdStats(.self$IdentPeptideData), silent=TRUE)
                        if (inherits(identStats, "try-error")) {
                            warning("No Regular and Random protein database types in identification peptide data.")
                        } else {
                            .self$IdentPeptideData <- cbind(.self$IdentPeptideData, identStats)
                            .self$SynapterLog <- c(.self$SynapterLog,
                                                   "Added identification statistics to identification data")
                            .self$filterIdentRandomEntries()
                        }
                    },
                    mergePeptides = function () {
                        'Merge quantitation and identification final peptide data'
                        .self$MergedFeatures <- merge(.self$IdentPeptideData,
                                                      .self$QuantPeptideData,
                                                      by.x = "peptide.seq",
                                                      by.y = "peptide.seq",
                                                      suffixes = c(".ident", ".quant"))
                        .self$MergedFeatures$deltaRt <-
                            .self$MergedFeatures$precursor.retT.ident -
                                .self$MergedFeatures$precursor.retT.quant
                        if (all(.self$MergedFeatures$deltaRt == 0)) {
                            stop("Merged identification and quantitation data have identical retention times. Modelling not possible")
                        }
                        if (any(.self$MergedFeatures$peptide.mhp.quant !=
                                .self$MergedFeatures$peptide.mhp.ident)) {
                            .self$MergedFeatures <- data.frame()
                            stop("Identification and quantitation theoretical peptide masses do not match.")
                        }
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste0("Merged identification and quantitation data [",
                                                      paste(dim(.self$MergedFeatures), collapse=","), "]"))
                    },
                    modelRetentionTime = function(span) {
                        'Models retention time'

                        if (missing(span)) {
                            span <- .self$LowessSpan
                        } else {
                           .self$LowessSpan <- span
                        }
                        .self$RtModel <- modelRetTime(.self$MergedFeatures$precursor.retT.ident,
                                                      .self$MergedFeatures$deltaRt,
                                                      span = span)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Modelled retention time using lowess and span ",
                                                     .self$LowessSpan, sep=""))
                        ## new 0.4.6 - use model to ad precited rt and sd
                        ##             to original IdentPeptideData
                        newIdentData <- doHDMSePredictions(.self$IdentPeptideData,
                                                           .self$RtModel) ## nsd missing
                        .self$IdentPeptideData$predictedRt <- newIdentData$predicted
                        .self$IdentPeptideData$sdRt <- newIdentData$sd
                    },
                    findEMRTs = function() {
                        if (length(.self$RtModel) == 0)
                            stop("First build a retention time model using 'modelRt'.")
                        if (length(.self$PpmError) == 0) {
                            warning("Mass error tolerance for EMRT matching undefined. Setting to default value.")
                            .self$setPpmError()
                        }
                        if (length(.self$RtNsd) == 0) {
                            warning("Retention time tolerance for EMRT matching undefined. Setting to default value.")
                            .self$setRtNsd()
                        }
                        if (length(.self$ImDiff) == 0) {
                            warning("Ion mobility difference for EMRT matching undefined. Setting to default value.")
                            .self$setImDiff()
                        }

                        .self$MatchedEMRTs <- findMSeEMRTs(.self$IdentPeptideData,
                                                           .self$QuantPep3DData,
                                                           .self$MergedFeatures,
                                                           .self$RtNsd,
                                                           .self$PpmError,
                                                           .self$ImDiff,
                                                           .self$RtModel)
                        .self$MatchedEMRTs$synapterPlgsAgreement <-
                          .findSynapterPlgsAgreement(.self$MatchedEMRTs)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste0("Matched identification peptides and quantitation EMRTs [",
                                                     paste(dim(.self$MatchedEMRTs), collapse = ","),
                                                     "]"))
                    },
                    rescueEMRTs = function(method) {
                        .self$MatchedEMRTs <- .rescueEMRTs(.self$MatchedEMRTs,
                                                           .self$MergedFeatures,
                                                           .self$QuantPep3DData,
                                                           method)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste0("Rescue EMRTs [",
                                                     paste(dim(.self$MatchedEMRTs), collapse = ","),
                                                     "] (", method, ")"))
                    },
                    fragmentMatching = function(verbose=TRUE) {
                      if (!nrow(.self$MatchedEMRTs)) {
                        stop("You have to run ", sQuote("findEMRTs"),
                             " first!")
                      }
                      if ("FragmentMatching" %in% colnames(.self$MatchedEMRTs)) {
                        ## remove previous FragmentMatching results
                        .self$MatchedEMRTs$FragmentMatching <-
                        .self$MatchedEMRTs$FragmentMatchingDiff <-
                        .self$MatchedEMRTs$FragmentMatchingRank <- NULL
                      }
                      emrts <- flatMatchedEMRTs(.self$MatchedEMRTs,
                                                .self$QuantPep3DData,
                                                na.rm=FALSE,
                                                verbose=verbose)
                      if (!length(.self$FragmentMatchingPpmTolerance)) {
                        warning("FragmentMatching ppm tolerance undefined. ",
                                "Setting to default value.")
                        .self$setFragmentMatchingPpmTolerance()
                      }

                      .self$FragmentMatching <-
                        .fragmentMatching(flatEmrts=emrts,
                                          identFragments=.self$IdentFragmentData,
                                          quantSpectra=.self$QuantSpectrumData,
                                          tolerance=.self$FragmentMatchingPpmTolerance/1e6,
                                          verbose=verbose)
                      .self$SynapterLog <-
                        c(.self$SynapterLog,
                          paste0("FragmentMatching of identification fragments ",
                                 "and quantitation spectra data at ",
                                 .self$FragmentMatchingPpmTolerance,
                                 "ppm [",
                                 paste0(dim(.self$FragmentMatching), collapse=","),
                                 "]"))
                      .self$MatchedEMRTs <- .appendFragmentMatchingColumn(
                        .self$MatchedEMRTs, .self$FragmentMatching)

                      .self$SynapterLog <-
                        c(.self$SynapterLog,
                          paste0("Append FragmentMatching results to ",
                                 "MatchedEMRTs [",
                                 paste0(dim(.self$MatchedEMRTs), collapse=","),
                                 "]"))
                    }))


## Synapter$lock("Version") ## this blocks $copy

## Grid search
.Synapter$methods(
                  list(
                       searchGrid = function(ppms, nsds, imdiffs, subset, n,
                                             verbose = TRUE) {
                         'Performs a grid search in ppm x nsd x imdiff space.'
                         .IdentPeptideData <- .self$IdentPeptideData
                         if (!missing(subset)) {
                           if (subset < 1) {
                             ns <- ceiling(nrow(.IdentPeptideData) * subset)
                             hidx <- sample(nrow(.IdentPeptideData), ns)
                             .IdentPeptideData <- .IdentPeptideData[hidx, ]
                             midx <- which(.self$MergedFeatures$precursor.leID.ident %in%
                                           .IdentPeptideData$precursor.leID)
                             if (length(hidx) < length(midx))
                               stop("Less identification peptides that merged ones!")
                           } else {
                             midx <- TRUE ## subset == 1
                           }
                           .log <- paste0("Performed grid search using ",
                                          subset * 100, "% of identification peptides.")
                         } else if (!missing(n)){
                           if (n > nrow(.IdentPeptideData))
                             warning("There are less than ", n, " peptides available. ",
                                     "Using all ", nrow(.HDMEsPeptideData), " peptides.")
                           hidx <- sample(nrow(.IdentPeptideData), n)
                           .IdentPeptideData <- .IdentPeptideData[hidx, ]
                           midx <- which(.self$MergedFeatures$precursor.leID.ident %in%
                                         .IdentPeptideData$precursor.leID)
                           .log <- paste0("Performed grid search using ",
                                          n, " identification peptides.")
                         }

                         ## test for ion mobility
                         if (!"precursor.Mobility" %in% colnames(.self$IdentPeptideData) ||
                             !"clust_drift" %in% colnames(.self$QuantPep3DData)) {
                           if (verbose) {
                             message("No ion mobility available. ",
                                     "Step back to 2D grid search.")
                           }
                           ## by setting imdiffs to Inf we disable the 3D grid
                           ## search
                           imdiffs <- Inf
                         }

                         .grid <- gridSearch3(.self$RtModel,
                                              .IdentPeptideData,
                                              .self$QuantPep3DData,
                                              .self$MergedFeatures[midx, ],
                                              ppms,
                                              nsds,
                                              imdiffs,
                                              verbose = verbose)
                         .self$Grid <- .grid[1:2]
                         .self$GridDetails <- .grid[[3]]
                         ## generate a proper detail grid - new v 0.7.2
                         .nd <- dim(.grid[[1]])
                         .dd <- do.call(rbind, .grid[[3]])
                         .dd[is.na(.dd)] <- 0
                         .prec <- .dd[,"1"]/(.dd[,"-1"]+.dd[,"1"])
                         .self$Grid[[3]] <- aperm(array(.prec,
                                                        dim=c(
                                                          length(imdiffs),
                                                          length(ppms),
                                                          length(nsds)),
                                                        dimnames=list(
                                                          as.character(imdiffs),
                                                          as.character(ppms),
                                                          as.character(nsds))),
                                                  perm=c(3, 2, 1))
                         names(.self$Grid)[3] <- "details"
                         .self$SynapterLog <- c(.self$SynapterLog, .log)
                       },
                       getBestGridValue = function() {
                         'Retieves the highest grid value.'
                         if ( length(.self$Grid) == 0 )
                           stop("No grid search result found.")
                         ans <- sapply(.self$Grid, max)
                         return(ans)
                       },
                       getBestGridParams = function() {
                         'Retrieves the best grid search (ppm, nsd, imdiff) triple(s).'
                         if ( length(.self$Grid) == 0 )
                           stop("No grid search result found.")
                         .getBestParams <- function(x) {
                           k <- arrayInd(which( x == max(x) ),
                                         dim(x),
                                         useNames = TRUE)
                           .dn <- dimnames(x)
                           k <- apply(k, 1,
                                      function(i) {
                                      as.numeric(c(.dn[[1]][i[1]],
                                                   .dn[[2]][i[2]],
                                                   .dn[[3]][i[3]]))
                                      })
                           rownames(k) <- c("nsd", "ppm", "imdiff")
                           return(t(k))
                         }
                         ans <- lapply(.self$Grid, .getBestParams)
                         return(ans)
                       },
                       setBestGridParams = function(what) {
                         'Sets best grid search (ppm, nsd, imdiff) triple.'
                         i <- 1 ## take first parameter pair by default
                         if (what == "auto") {
                           x <- .self$getBestGridParams()[["prcntModel"]]
                           if (nrow(x) > 1) {
                             y <- .self$getBestGridParams()[["details"]]
                             ppm_i <- x[, "ppm"] %in% y[, "ppm"]
                             nsd_i <- x[, "nsd"] %in% y[, "nsd"]
                             imdiff_i <- x[, "imdiff"] %in% y[, "imdiff"]
                             if (any(ppm_i & nsd_i & imdiff_i)) {
                               ## (1) first ppm and nsd match
                               i <- which(ppm_i & nsd_i & imdiff_i)[1]
                             } else if (any(ppm_i)) {
                               ## (2) first ppm match
                               i <- which(ppm_i)[1]
                             } else if (any(nsd_i)) {
                               ## (3) first nsd match
                               i <- which(nsd_i)[1]
                             } else if (any(imdiff_i)) {
                               ## (4) first imdiff match
                               i <- which(imdiff_i)[1]
                             }
                             ## else (5) no match - taking first one
                           }
                         } else {
                           x <- switch(what,
                                       total   = .self$getBestGridParams()[["prcntTotal"]],
                                       model   = .self$getBestGridParams()[["prcntModel"]],
                                       details = .self$getBestGridParams()[["details"]])
                         }
                         .self$RtNsd <- x[i,"nsd"]
                         .self$PpmError <- x[i,"ppm"]
                         .self$ImDiff <- x[i,"imdiff"]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste0("Set 'nsd', 'ppm error' and 'im diff' to ",
                                                       .self$RtNsd, ", ", .self$PpmError, " and ",
                                                       .self$ImDiff, " (best '", what, "')"))
                       }))

## Setters
.Synapter$methods(list(
                       setPepScoreFdr = function(fdr = 0.01) {
                         'Sets peptide score fdr.'
                         .self$PepScoreFdr <- fdr
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set peptide score fdr to ",
                                                      .self$PepScoreFdr,
                                                      sep=""))
                       },
                       setProtFpr = function(fpr = 0.01) {
                         'Sets protein fpr threshold.'
                         .self$ProtFpr <- fpr
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set protein fpr to ",
                                                      .self$ProtFpr,
                                                      sep=""))
                       },
                       setIdentPpmError = function(ppm = 10) {
                         'Sets identification mass error ppm threshold.'
                         .self$IdentPpmError <- ppm
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set identification ppm error to ",
                                                      .self$IdentPpmError,
                                                      sep=""))
                       },
                       setQuantPpmError = function(ppm = 10) {
                         'Sets quantitation mass error ppm threshold.'
                         .self$QuantPpmError <- ppm
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set quantitation ppm error to ",
                                                      .self$QuantPpmError,
                                                      sep=""))
                       },
                       setPpmError = function(ppm = 10) {
                         'Sets matching mass error tolerance threshold.'
                         .self$PpmError <- ppm
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set ppm error to ",
                                                      .self$PpmError,
                                                      sep=""))
                       },
                       setLowessSpan = function(span = 0.05) {
                         'Sets lowess\' span parameter.'
                         .self$LowessSpan <- span
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set span to ",
                                                      .self$LowessSpan,
                                                      sep=""))
                       },
                       setRtNsd = function(nsd = 2) {
                         .self$RtNsd <- nsd
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set nsd to ",
                                                      .self$RtNsd,
                                                      sep=""))
                       },
                       setImDiff = function(imdiff = 0.5) {
                         .self$ImDiff <- imdiff
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste0("Set imdiff to ",
                                                       .self$ImDiff))
                       },
                       setFragmentMatchingPpmTolerance = function(ppm = 25) {
                         'Sets FragmentMatching mass error tolerance threshold.'
                         .self$FragmentMatchingPpmTolerance <- ppm
                         .self$SynapterLog <-
                           c(.self$SynapterLog,
                             paste0("Set FragmentMatching ppm error to ", ppm))
                       }))

## GLOBAL FILTERS
.Synapter$methods(list(filterUniqueSeq = function() {
                         'Keep unique (non duplicated) identification and quantitation peptide sequences.'
                         ## Since v 0.5.0, filtering .self$[HD]QuantPeptideData dataframes (rather than mf ones)
                         ## when data is loaded, right after filtering for db unique peptides, to get rid of
                         ## pre/post-fix peptides (like ABCD and ABCDXXXXX), that are not filtered out with
                         ## filterUniqueDbPeptides
                         .self$filterUniqueQuantSeq()
                         .self$filterUniqueIdentSeq()
                       },
                       filterPeptideLength = function(l) {
                           'Filters quantitation and identificatio run peptides of length < l out.'
                           .filterPeptideLength <- function(seq, l)
                               nchar(seq) >= l
                           selQt <- .filterPeptideLength(.self$QuantPeptideData$peptide.seq, l)
                           selId <- .filterPeptideLength(.self$IdentPeptideData$peptide.seq, l)
                           .self$QuantPeptideData <- .self$QuantPeptideData[selQt, ]
                           .self$IdentPeptideData   <- .self$IdentPeptideData[selId, ]
                           .self$SynapterLog <- c(.self$SynapterLog,
                                                  paste("Keeping quant peptides of length >= ", l,
                                                        " [", paste(dim(.self$QuantPeptideData), collapse=","),
                                                        "]", sep = ""))
                           .self$SynapterLog <- c(.self$SynapterLog,
                                                  paste("Keeping ident peptides of length >= ", l,
                                                        " [", paste(dim(.self$IdentPeptideData), collapse=","),
                                                        "]", sep = ""))
                       },
                       filterUniqueQuantSeq = function() {
                         'Keep unique (non duplicated) quantitation peptide sequences.'
                         quant2keep <- names(which(table(.self$QuantPeptideData$peptide.seq) == 1))
                         .self$QuantPeptideData <-
                           .self$QuantPeptideData[.self$QuantPeptideData$peptide.seq %in% quant2keep, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Kept quantitation non duplicated peptide sequences [",
                                                      paste(dim(.self$QuantPeptideData),
                                                            collapse=","), "]", sep = ""))
                       },
                       filterUniqueIdentSeq = function() {
                         'Keep unique (non duplicated) identification peptide sequences.'
                         ident2keep <- names(which(table(.self$IdentPeptideData$peptide.seq) == 1))
                         .self$IdentPeptideData <-
                           .self$IdentPeptideData[.self$IdentPeptideData$peptide.seq %in% ident2keep, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Kept identification non duplicated peptide sequences [",
                                                      paste(dim(.self$IdentPeptideData),
                                                            collapse=","), "]", sep = ""))
                       },
                       filterQuantFunction = function() {
                         'Filters quantitation Pep3D data by Function == 1.'
                         .self$QuantPep3DData <- filterFunction(.self$QuantPep3DData)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered quantitation Pep3D by Function [",
                                                      paste(dim(.self$QuantPep3DData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterDuplicatedQuantSpectrumIds = function() {
                         'Removes duplicated quantitation EMRT spectrum ids (different charge states, isotopes,.. ) keeping the one with isFid == 1 (is used for identification).'
                         keep  <- .self$QuantPep3DData$isFid == 1
                         .self$QuantPep3DData <- .self$QuantPep3DData[keep, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Kept unique spectrum ids ",
                                                      "quantitation Pep3D [", paste(dim(.self$QuantPep3DData), collapse=","), "] ",
                                                      sep=""))
                       },
                       filterQuantMatchType = function() {
                         'Filter quantitation peptide file for PepFrag1 and PepFrag2'
                         .self$QuantPeptideData <- filterPeptideMatchType(.self$QuantPeptideData)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered quantitation peptide by PepFrag1|2 [",
                                                      paste(dim(.self$QuantPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterIdentMatchType = function() {
                         'Filter identification peptide file for PepFrag1 and PepFrag2'
                         .self$IdentPeptideData <- filterPeptideMatchType(.self$IdentPeptideData)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered identification peptide by PepFrag1|2 [",
                                                      paste(dim(.self$IdentPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterMatchType = function() {
                         'Filter identification and quantitation peptide file for PepFrag1 and PepFrag2 peptides.'
                         .self$filterIdentMatchType()
                         .self$filterQuantMatchType()
                       },
                       filterIdentRandomEntries = function() {
                         .self$IdentPeptideData <-
                           .self$IdentPeptideData[.self$IdentPeptideData$protein.dataBaseType == "Regular", ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered identification Random entries [",
                                                      paste(dim(.self$IdentPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterQuantRandomEntries = function() {
                         .self$QuantPeptideData <-
                           .self$QuantPeptideData[.self$QuantPeptideData$protein.dataBaseType == "Regular", ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered quantitation Random entries [",
                                                      paste(dim(.self$QuantPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterRandomEntries = function() {
                         .self$filterIdentRandomEntries()
                         .self$filterQuantRandomEntries()
                       },
                       filterMismatchingQuantIntensities = function() {
                        ## this method should not be necessary; a mismatch
                        ## between intensity values should never happen.
                        ## See https://github.com/lgatto/synapter/issues/42
                        ## for details
                        idx <- match(.self$QuantPeptideData$precursor.leID,
                                     .self$QuantPep3DData$spectrumID)

                        iMismatch <- which(.self$QuantPeptideData$precursor.inten !=
                                           .self$QuantPep3DData$Counts[idx])
                        nMismatch <- length(iMismatch)

                        if (nMismatch == length(idx)) {
                          stop("It seems that all IDs are correct but ",
                               "there are not any matching intensity values.")
                        }
                        if (nMismatch) {
                          warning("Filtering ", nMismatch, " (of ", length(idx), " total) ",
                                  "entries of the quantitation final peptide and Pep3D file because ",
                                  "they differ in their intensity values.")
                          .self$QuantPeptideData <- .self$QuantPeptideData[-iMismatch, ]
                          .self$QuantPep3DData <- .self$QuantPep3DData[-idx[iMismatch], ]
                          .self$SynapterLog <- c(.self$SynapterLog,
                                                 paste0("Filtered ", nMismatch,
                                                        " quantitation final peptide data entries because of intensity mismatch [",
                                                        paste(dim(.self$QuantPeptideData), collapse = ","), "]"))
                          .self$SynapterLog <- c(.self$SynapterLog,
                                                 paste0("Filtered ", nMismatch,
                                                        " quantitation Pep3D data entries because of intensity mismatch [",
                                                       paste(dim(.self$QuantPep3DData), collapse = ","), "]"))
                        }
                       }
                       ))


## PEPTIDE DATA FILTER
.Synapter$methods(list(
                       filterQuantPepScore = function(fdrMethod, log = TRUE) {
                         'Filter quantitation peptide data based on peptide score fdr.'
                         if (length(.self$PepScoreFdr) == 0) {
                           warning("Pep score FDR undefined. Setting to default value.")
                           .self$setPepScoreFdr()
                         }
                         .self$QuantPeptideData <- filterPepScore(.self$QuantPeptideData,
                                                                .self$PepScoreFdr,
                                                                fdrMethod)
                         if (log)
                           .self$SynapterLog <- c(.self$SynapterLog,
                                                  paste0("Filtered quantitation peptide data on pep score with fdr ",
                                                        .self$PepScoreFdr, " [",
                                                        paste(dim(.self$QuantPeptideData), collapse=","), "]"))
                       },
                       filterIdentPepScore = function(fdrMethod, log = TRUE) {
                         'Filter identification peptide data based on peptide score fdr.'
                         if (length(.self$PepScoreFdr) == 0) {
                           warning("Pep score FDR undefined. Setting to default value.")
                           .self$setPepScoreFdr()
                         }
                         .self$IdentPeptideData <- filterPepScore(.self$IdentPeptideData,
                                                                  .self$PepScoreFdr,
                                                                  fdrMethod)
                         if (log)
                           .self$SynapterLog <- c(.self$SynapterLog,
                                                  paste("Filtered identification peptide data on pep score with fdr ",
                                                        .self$PepScoreFdr, " [",
                                                        paste(dim(.self$IdentPeptideData), collapse=","),
                                                        "]", sep=""))
                       },
                       filterIdentPpmError = function(ppm = .self$IdentPpmError) {
                         'Filter identification mass error.'
                         if (length(ppm) == 0) {
                           warning("Identification ppm error undefined. Setting to default value.")
                           .self$setIdentPpmError()
                           ppm <- .self$IdentPpmError
                         }
                         .self$IdentPeptideData <-
                           filter.error.ppm(.self$IdentPeptideData,
                                            colname = "errorppm",
                                            ppm = ppm)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Identification mass error filtered at ", .self$IdentPpmError, " ppm ",
                                                      "[", paste(dim(.self$IdentPeptideData), collapse=","), "]",
                                                      sep=""))
                       },
                       filterQuantPpmError = function(ppm = .self$QuantPpmError) {
                         'Filter quantitation mass error.'
                         if (length(ppm) == 0) {
                           warning("quantitation ppm error undefined. Setting to default value.")
                           .self$setQuantPpmError()
                           ppm <- .self$QuantPpmError
                         }
                         .self$QuantPeptideData <-
                           filter.error.ppm(.self$QuantPeptideData,
                                            colname = "errorppm",
                                            ppm = ppm)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Quantitation mass error filtered at ", .self$QuantPpmError, " ppm ",
                                                      "[", paste(dim(.self$QuantPeptideData), collapse=","), "]",
                                                      sep=""))
                       },
                       filterIdentProtFpr = function() {
                         'Filter identification peptide data based on protein fpr.'
                         if (length(.self$ProtFpr) == 0) {
                           warning("Protein FPR undefined. Setting to default value.")
                           .self$setProtFpr()
                         }
                         .self$IdentPeptideData <- filterProtFpr(.self$IdentPeptideData,
                                                                   .self$ProtFpr)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered identification peptide data using protein fpr ",
                                                      .self$ProtFpr, " [",
                                                      paste(dim(.self$IdentPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       filterQuantProtFpr = function() {
                         'Filter quantitation peptide data based on protein fpr.'
                         if (length(.self$ProtFpr) == 0) {
                           warning("Protein FPR undefined. Setting to default value.")
                           .self$setProtFpr()
                         }
                         .self$QuantPeptideData <- filterProtFpr(.self$QuantPeptideData,
                                                                 .self$ProtFpr)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Filtered quantitation peptide data using protein fpr ",
                                                      .self$ProtFpr, " [",
                                                      paste(dim(.self$QuantPeptideData), collapse=","),
                                                      "]", sep=""))
                       },

                       filterUniqueQuantDbPeptides = function(filename, missedCleavages = 0,  IisL = FALSE, verbose = TRUE) {
                         'Filters quantitation tryptic peptides that match one and only one protein in the fasta database.'
                         .self$filterUniqueDbPeptides(filename,
                                                      what="quant",
                                                      missedCleavages=missedCleavages,
                                                      IisL=IisL,
                                                      verbose=verbose)
                       },

                       filterUniqueIdentDbPeptides = function(filename, missedCleavages = 0, IisL = FALSE, verbose = TRUE) {
                         'Filters identification tryptic peptides that match one and only one protein in the fasta database.'
                         .self$filterUniqueDbPeptides(filename,
                                                      what="ident",
                                                      missedCleavages=missedCleavages,
                                                      IisL=IisL,
                                                      verbose=verbose)
                       },

                       filterUniqueDbPeptides = function(filename, what = c("ident", "quant"), missedCleavages = 0, IisL = FALSE, verbose = TRUE) {
                        'Filters tryptic peptides that match one and only one protein in the fasta database.'
                        what <- match.arg(what, several.ok=TRUE)

                        upepset <- dbUniquePeptideSet(filename,
                                                      missedCleavages=missedCleavages,
                                                      IisL=IisL,
                                                      verbose=verbose)
                        .self$DbFastaFile <- filename
                        logmsg <- paste0("peptides that match unique protein (",
                                         attr(upepset, "missedCleavages"),
                                         " missed cleavages; I/L treatment: ",
                                         ifelse(attr(upepset, "IisL"), "I == L", "I != L"), ")")

                        if ("ident" %in% what) {
                          sel <- .self$IdentPeptideData$peptide.seq %in% upepset
                          .self$IdentPeptideData <- .self$IdentPeptideData[sel, ]
                          .self$SynapterLog <- c(.self$SynapterLog,
                                                 paste0("Kept identification ", logmsg,
                                                       "[", paste(dim(.self$IdentPeptideData), collapse=","), "]"))
                        }

                        if ("quant" %in% what) {
                          sel <- .self$QuantPeptideData$peptide.seq %in% upepset
                         .self$QuantPeptideData <- .self$QuantPeptideData[sel, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste0("Kept quantitation ", logmsg,
                                                       "[", paste(dim(.self$QuantPeptideData), collapse=","), "]"))
                        }
                       },

                       filterFragments = function(what, minIntensity = NULL,
                                                  maxNumber = NULL, verbose = TRUE) {
                         'Filters spectra/fragments using a minimal intensity or a maximal number of fragments as threshold.'

                         what <- match.arg(what, choices = c("fragments.ident",
                                                             "spectra.quant"))
                         msexp <- switch(what,
                            "fragments.ident" = .self$IdentFragmentData,
                            "spectra.quant" = .self$QuantSpectrumData)

                         if (!length(msexp)) {
                           stop("You have to import the ", sQuote(what),
                                " data first!")
                         }

                         msg <- "Filtered "

                         if (what == "fragments.ident") {
                           .self$IdentFragmentData <-
                             .filterIntensity(msexp,
                                              minIntensity = minIntensity,
                                              maxNumber = maxNumber,
                                              verbose = verbose)
                           msg <- paste0("identification fragment data")
                         }

                         if (what == "spectra.quant") {
                           .self$QuantSpectrumData <-
                             .filterIntensity(msexp,
                                              minIntensity = minIntensity,
                                              maxNumber = maxNumber,
                                              verbose = verbose)
                           msg <- paste0("quantitation spectra")
                         }

                         msg <- paste0(" using a ",
                                       ifelse(is.null(minIntensity), "",
                                              paste0("minimal intensity > ",
                                                     minIntensity)),
                                       ifelse(is.null(maxNumber), "",
                                              paste0("maximal number < ",
                                                     maxNumber)))

                         .self$SynapterLog <- c(.self$SynapterLog, msg)
                       },

                       filterUniqueMatches = function(minNumber) {
                         'Filters unique matches using FragmentMatching results.'

                         if (!nrow(.self$FragmentMatching)) {
                           stop("You have to run ", sQuote("fragmentMatching"),
                                " first!")
                         }

                         .self$MatchedEMRTs <-
                           .filterUniqueMatches(obj = .self,
                                                mincommon = minNumber)
                         .self$MatchedEMRTs$synapterPlgsAgreement <-
                           .findSynapterPlgsAgreement(.self$MatchedEMRTs)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste0("Filtered unique matched EMRTs ",
                                                       "(minimal number of common peaks: ",
                                                       minNumber, " [",
                                                      paste0(dim(.self$MatchedEMRTs),
                                                             collapse=","), "]"))
                       },

                       filterNonUniqueMatches = function(minDelta) {
                         'Filters non unique matches using FragmentMatching results.'

                         if (!nrow(.self$FragmentMatching)) {
                           stop("You have to run ", sQuote("fragmentMatching"),
                                " first!")
                         }

                         .self$MatchedEMRTs <-
                           .filterNonUniqueMatches(obj = .self,
                                                   mindelta = minDelta)
                         .self$MatchedEMRTs$synapterPlgsAgreement <-
                           .findSynapterPlgsAgreement(.self$MatchedEMRTs)
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste0("Filtered non unique matched EMRTs ",
                                                       "(minimal delta of common peaks: ",
                                                       minDelta, " [",
                                                      paste0(dim(.self$MatchedEMRTs),
                                                             collapse=","), "]"))
                       }))
