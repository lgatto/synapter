.Synapter <-
    setRefClass("Synapter",
                fields = list(
                    ## introspection
                    Version = "character",
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
                    MatchedEMRTs = "data.frame"),
                methods = list(
                    initialize = function() {
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
                    },
                    loadMasterData = function() {
                        message("Reading quantitation final peptide file...")
                        .self$QuantPeptideData <- read.csv(.self$QuantPeptideFile, stringsAsFactors = FALSE)
                        .self$QuantPeptideData$errorppm <-
                            error.ppm(obs = .self$QuantPeptideData$precursor.mhp,
                                      theo = .self$QuantPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation peptide data [",
                                                     paste(dim(.self$QuantPeptideData), collapse = ","),
                                                     "]", sep=""))
                        if (!.self$Master)
                            stop("Identification final peptide is not a master file")
                        message("Reading master identification peptide file...")
                        ext <- getExtension(.self$IdentPeptideFile)
                        if (tolower(ext) == "csv") {
                            .self$IdentPeptideData <- read.csv(.self$IdentPeptideFile,
                                                               stringsAsFactors = FALSE)
                            .self$SynapterLog <- c(.self$SynapterLog,
                                                   paste("Read master identification peptide (csv) data [",
                                                         paste(dim(.self$IdentPeptideData), collapse = ","),
                                                         "]", sep=""))
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
                        } else {
                            stop("Master peptide file extention not recognised. Must be 'csv' or 'rds'.")
                        }
                        message("Reading quantitation Pep3D file...")
                        .self$QuantPep3DData <- read.csv(.self$QuantPep3DFile, stringsAsFactors = FALSE)
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
                        .self$filterDuplicatedQuantSpectrumIds()
                        message("Computing quantitation identification statistics...")
                        .self$addQuantIdStats() ## Quant only
                    },                           
                    loadData = function() {
                        if (.self$Master)
                            stop("Identification final peptide is a master file")
                        message("Reading identification final peptide file...")
                        .self$IdentPeptideData <- read.csv(.self$IdentPeptideFile, stringsAsFactors = FALSE)
                        .self$IdentPeptideData$errorppm <-
                            error.ppm(obs = .self$IdentPeptideData$precursor.mhp,
                                      theo = .self$IdentPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read identification peptide data [",
                                                     paste(dim(.self$IdentPeptideData), collapse = ","),
                                                     "]", sep=""))
                        message("Reading quantitation final peptide file...")
                        .self$QuantPeptideData <- read.csv(.self$QuantPeptideFile, stringsAsFactors = FALSE)
                        .self$QuantPeptideData$errorppm <-
                            error.ppm(obs = .self$QuantPeptideData$precursor.mhp,
                                      theo = .self$QuantPeptideData$peptide.mhp)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Read quantitation peptide data [",
                                                     paste(dim(.self$QuantPeptideData), collapse = ","),
                                                     "]", sep=""))
                        message("Reading quantitation Pep3D file...")
                        .self$QuantPep3DData <- read.csv(.self$QuantPep3DFile, stringsAsFactors = FALSE)
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
                        .self$filterDuplicatedQuantSpectrumIds()
                        message("Computing identification statistics...")
                        .self$addIdStats()
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
                        if (missing(span))
                            span <- .self$LowessSpan
                        .self$RtModel <- modelRetTime(.self$MergedFeatures, span = span)
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
                    findEMRTs = function(mergedEMRTs) {
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
                        .self$MatchedEMRTs <- findMSeEMRTs(.self$IdentPeptideData,
                                                           .self$QuantPep3DData,
                                                           .self$MergedFeatures,
                                                           .self$RtNsd,
                                                           .self$PpmError,
                                                           .self$RtModel,
                                                           mergedEMRTs)
                        .self$SynapterLog <- c(.self$SynapterLog,
                                               paste("Matched identification peptides and quantitation EMRTs [",
                                                     paste(dim(.self$MatchedEMRTs), collapse = ","),
                                                     "] (", mergedEMRTs, ")",
                                                     sep = ""))
                    }))
                  

## Synapter$lock("Version") ## this blocks $copy

## Grid search
.Synapter$methods(
                  list(
                       searchGrid = function(ppms, nsds, subset, n, verbose = TRUE) {
                         'Performs a grid search in ppm x nsd space.'
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
                         .grid <- gridSearch2(.self$RtModel,
                                              .IdentPeptideData,
                                              .self$QuantPep3DData,
                                              .self$MergedFeatures[midx, ],
                                              ppms,
                                              nsds,
                                              verbose = verbose)
                         .self$Grid <- .grid[1:2]
                         .self$GridDetails <- .grid[[3]]
                         ## generate a proper detail grid - new v 0.7.2
                         .nc <- ncol(.grid[[1]])
                         .nr <- nrow(.grid[[1]])
                         .gdet <- matrix(sapply(.self$GridDetails,
                                                function(.x) {
                                                  .x1 <- ifelse(is.na(.x["1"]) , 0, .x["1"])
                                                  .x2 <- ifelse(is.na(.x["-1"]), 0, .x["-1"])
                                                  ans <- .x1/(.x1 + .x2)
                                                  if (is.na(ans))
                                                    ans <- 0
                                                  return(ans)
                                                }),
                                         ncol = .nc,
                                         nrow = .nr,
                                         byrow = TRUE)
                         rownames(.gdet) <- rownames(.grid[[1]])
                         colnames(.gdet) <- colnames(.grid[[1]])
                         .self$Grid[[3]] <- .gdet
                         names(.self$Grid)[3] <- "details"
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                .log)
                       },
                       getBestGridValue = function() {
                         'Retieves the highest grid value.'
                         if ( length(.self$Grid) == 0 )
                           stop("No grid search result found.")
                         ans <- sapply(.self$Grid, max)
                         return(ans)
                       },
                       getBestGridParams = function() {
                         'Retrieves the best grid search (ppm, nsd) pair(s).'
                         if ( length(.self$Grid) == 0 )
                           stop("No grid search result found.")
                         .getBestParams <- function(x) {
                           k <- arrayInd(which( x == max(x) ),
                                         dim(x),
                                         useNames = TRUE)
                           k <- apply(k, 1,
                                      function(i)
                                      as.numeric(c(colnames(x)[i["col"]], rownames(x)[i["row"]])))
                           rownames(k) <- c("ppm", "nsd")
                           return(t(k))
                         }
                         ans <- lapply(.self$Grid, .getBestParams)
                         return(ans)
                       },
                       setBestGridParams = function(what) {
                         'Sets best grid search (ppm, nsd) pair.'
                         i <- 1 ## take first parameter pair by default
                         if (what == "auto") {
                           x <- .self$getBestGridParams()[["prcntModel"]]
                           if (nrow(x) > 1) {
                             y <- .self$getBestGridParams()[["details"]]
                             ppm_i <- x[, "ppm"] %in% y[, "ppm"]
                             nsd_i <- x[, "nsd"] %in% y[, "nsd"]
                             if (any(ppm_i & nsd_i)) {
                               ## (1) first ppm and nsd match
                               i <- which(ppm_i & nsd_i)[1]
                             } else if (any(ppm_i)) {
                               ## (2) first ppm match
                               i <- which(ppm_i)[1]
                             } else if (any(ppm_i)) {
                               ## (2) first nsd match
                               i <- which(nsd_i)[1]
                             } 
                             ## else (4) no match - taking first one
                           }
                         } else {
                           x <- switch(what,
                                       total   = .self$getBestGridParams()[["prcntTotal"]],
                                       model   = .self$getBestGridParams()[["prcntModel"]],
                                       details = .self$getBestGridParams()[["details"]])
                         }                         
                         .self$RtNsd <- x[i,"nsd"]
                         .self$PpmError <- x[i,"ppm"]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Set 'nsd' and 'ppm error' to ",
                                                      .self$RtNsd, " and ", .self$PpmError,
                                                      " (best '", what, "')", sep=""))
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
                         'Removes duplicated quantitation EMRT spectrum ids (different charge states, isotopes,.. ) keeping the first instance.'
                         keep <- !duplicated(.self$QuantPep3DData$spectrumID)
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

                       filterUniqueQuantDbPeptides = function(filename, missedCleavages = 0, verbose = TRUE) {
                         'Filters quantitation tryptic peptides that match one and only one protein in the fasta database.'
                         upepset <- dbUniquePeptideSet(filename, missedCleavages, verbose)
                         .self$DbFastaFile <- filename
                         sel2 <- .self$QuantPeptideData$peptide.seq %in% upepset
                         .self$QuantPeptideData <- .self$QuantPeptideData[sel2, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Kept quantitation peptides that match unique protein (",
                                                      missedCleavages, " missed cleavages) " ,
                                                      "[", paste(dim(.self$QuantPeptideData), collapse=","),
                                                      "]", sep=""))
                       },
                       
                       filterUniqueIdentDbPeptides = function(filename, missedCleavages, verbose = TRUE) {
                         'Filters identification tryptic peptides that match one and only one protein in the fasta database.'
                         upepset <- dbUniquePeptideSet(filename, missedCleavages, verbose)
                         .self$DbFastaFile <- filename
                         sel1 <- .self$IdentPeptideData$peptide.seq %in% upepset
                         .self$IdentPeptideData <- .self$IdentPeptideData[sel1, ]
                         .self$SynapterLog <- c(.self$SynapterLog,
                                                paste("Kept identification peptides that match unique protein (",
                                                      missedCleavages, " missed cleavages) " ,
                                                      "[", paste(dim(.self$IdentPeptideData), collapse=","),
                                                      "]", sep=""))
                       },

                       filterUniqueDbPeptides = function(filename, missedCleavages, verbose = TRUE) {
                         'Filters tryptic peptides that match one and only one protein in the fasta database.'
                         .self$filterUniqueIdentDbPeptides(filename, missedCleavages, verbose = verbose)
                         .self$filterUniqueQuantDbPeptides(filename, missedCleavages, verbose = FALSE)
                       }))
                       
