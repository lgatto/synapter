context("master")

test_that(".calculatePeptideFileCombinations", {
  r1 <- list(c(1, 2), c(1, 3), c(2, 3), 1:3)
  r2 <- list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4),
             c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4), 1:4)
  r3 <- list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))
  expect_equal(synapter:::.calculatePeptideFileCombinations(3), r1)
  expect_equal(synapter:::.calculatePeptideFileCombinations(4), r2)
  expect_equal(synapter:::.calculatePeptideFileCombinations(4, 2), r3)
  expect_silent(synapter:::.calculatePeptideFileCombinations(4, 2))
  expect_message(synapter:::.calculatePeptideFileCombinations(3, verbose=TRUE),
                 "3 peptide files available - 4 combinations")
  expect_message(synapter:::.calculatePeptideFileCombinations(4, verbose=TRUE),
                 "4 peptide files available - 11 combinations")
  expect_message(synapter:::.calculatePeptideFileCombinations(4, 2, verbose=TRUE),
                 "4 peptide files available - 6 combinations (limited by maxFileComb=2)*")
})

test_that(".peptideMatrix", {
  l <- list(LETTERS[1:5], LETTERS[1:3], LETTERS[3:7])
  m <- matrix(c(rep(1L, 5), 0L, 0L,
                rep(1L, 3), rep(0L, 4),
                0L, 0L, rep(1L, 5)), nrow = 7, ncol = 3,
              dimnames = list(LETTERS[1:7], 1:3))
  expect_equal(synapter:::.peptideMatrix(l), m)
})

test_that(".masterFdrSummary", {
  incorrect <- c(4, 5, 4)
  uni <- c(5, 7, 7)
  r <- data.frame(incorrect = incorrect, unique = uni,
                  proteotypic = uni,
                  fdr = incorrect / uni)
  r$combination <- list(c(1, 2), c(1, 3), c(2, 3))
  r$nbSample <- rep(2, 3)
  expect_equal(synapter:::.masterFdrSummary(r$combination, LETTERS[1:10],
                                            list(LETTERS[1:5],
                                                 LETTERS[1:3],
                                                 LETTERS[3:7]), fdr = 0.5), r)
})

test_that(".orderForMasterModels", {
  l <- list(A=list(IdentPeptideData=data.frame(A=1:2)),
            B=list(IdentPeptideData=data.frame(B=3:10)),
            C=list(IdentPeptideData=data.frame(C=3:5)))
  expect_equal(synapter:::.orderForMasterModels(l),
               list(c(2, 3, 1), c(3, 1, 2)))
})

test_that(".filterDuplicatedPeptideSequences", {
  l <- list(A=list(IdentPeptideData=data.frame(A=1:4,
                                        peptide.seq=LETTERS[rep(1:2, 2)],
                                        stringsAsFactors=FALSE)),
            B=list(IdentPeptideData=data.frame(A=1:6,
                                        peptide.seq=LETTERS[c(1, 1, 1, 2, 2, 3)],
                                        stringsAsFactors=FALSE)))
  r <- list(A=list(IdentPeptideData=data.frame(A=1:2,
                                        peptide.seq=LETTERS[1:2],
                                        stringsAsFactors=FALSE)),
            B=list(IdentPeptideData=data.frame(A=c(1, 4, 6),
                                        peptide.seq=LETTERS[1:3],
                                        stringsAsFactors=FALSE)))
  rownames(r$B$IdentPeptideData) <- c(1L, 4L, 6L)
  expect_equal(synapter:::.filterDuplicatedPeptideSequences(l[[1]]), r[[1]])
  expect_equal(synapter:::.filterDuplicatedPeptideSequences(l[[2]]), r[[2]])
})

test_that(".mergePeptideData", {
  l <- list(A=data.frame(peptide.seq = LETTERS[1:5],
                         precursor.retT = 1:5,
                         precursor.inten = 2*(1:5),
                         stringsAsFactors = FALSE),
            B=data.frame(peptide.seq = LETTERS[c(1:3, 6:10)],
                         precursor.retT = c(2:3, 10, 1:5),
                         precursor.inten = 1:8,
                         stringsAsFactors = FALSE))
  m <- data.frame(peptide.seq = LETTERS[1:3],
                  precursor.retT.master = 1:3,
                  precursor.inten.master = 2*(1:3),
                  precursor.retT.slave = c(2:3, 10),
                  precursor.inten.slave = 1:3,
                  deltaRt = c(-1, -1, -7),
                  intenRatio = 1,
                  stringsAsFactors = FALSE)

  expect_equal(synapter:::.mergePeptideData(l$A, l$B, verbose = FALSE), m)
  expect_equal(synapter:::.mergePeptideData(l$A, l$B, maxDeltaRt = 5,
                                            verbose = FALSE), m[-3,])
  expect_message(synapter:::.mergePeptideData(l$A, l$B, maxDeltaRt = 5,
                                              verbose = TRUE),
                 "Ignoring 1 features because deltaRt > 5")
})

test_that(".mergeMaster", {
  l <- list(A=list(IdentPeptideData=data.frame(
               peptide.seq = LETTERS[c(1:3, 6:7)],
               precursor.inten = 2*c(1:3, 6:7),
               precursor.retT = c(0.2, 0.4, 0.6, 1.2, 1.4),
               stringsAsFactors = FALSE)),
            B=list(IdentPeptideData=data.frame(
               peptide.seq = LETTERS[c(1:2, 4)],
               precursor.inten = c(1:2, 4),
               precursor.retT = c(0.21, 0.41, 0.81),
               stringsAsFactors = FALSE)),
            C=list(IdentPeptideData=data.frame(
               peptide.seq = LETTERS[c(1:2, 5)],
               precursor.inten = c(1:2, 5),
               precursor.retT = c(0.22, 0.42, 1.02),
               stringsAsFactors = FALSE)))
  m <- data.frame(peptide.seq = LETTERS[c(1:3, 6:7, 4:5)],
                  precursor.inten = c(2, 4, 6, 12, 14, 4, 5),
                  precursor.retT = c(0.2, 0.4, 0.6, 1.2, 1.4, 0.8, 1.0),
                  stringsAsFactors = FALSE)
  expect_error(synapter:::.mergeMaster(l[1]),
               "To create a master at least two identification data are needed.")
  expect_equal(suppressWarnings(
                 synapter:::.mergeMaster(l, span=0.5, verbose = FALSE)), m,
               tolerance = 0.02)
})

test_that(".regeneratePrecursorLeId", {
  l <- list(A=data.frame(peptide.seq = LETTERS[1:10],
                         precursor.leID = 10:1,
                         stringsAsFactors = FALSE),
            B=data.frame(peptide.seq = LETTERS[5:1],
                         precursor.leID = 1:5,
                         stringsAsFactors = FALSE))
  r <- list(A=data.frame(peptide.seq = LETTERS[1:10],
                         precursor.leID = 1:10,
                         stringsAsFactors = FALSE),
            B=data.frame(peptide.seq = LETTERS[5:1],
                         precursor.leID = 5:1,
                         stringsAsFactors = FALSE))
  expect_equal(synapter:::.regeneratePrecursorLeId(l, verbose = FALSE), r)
  expect_message(synapter:::.regeneratePrecursorLeId(l, verbose = TRUE),
                 "Regenerate precursor.leIDs")
})
