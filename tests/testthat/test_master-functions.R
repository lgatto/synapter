context("master")

test_that(".orderForMasterModels", {
    l <- list(A=list(IdentPeptideData=data.frame(A=1:2)),
              B=list(IdentPeptideData=data.frame(B=3:10)),
              C=list(IdentPeptideData=data.frame(C=3:5)))
    expect_equal(synapter:::.orderForMasterModels(l),
                 list(c(2, 3, 1), c(3, 1, 2)))
})

test_that(".mergePeptideData", {
    l <- list(A=data.frame(peptide.seq = LETTERS[1:5],
                           precursor.retT = 1:5,
                           stringsAsFactors = FALSE),
              B=data.frame(peptide.seq = LETTERS[c(1:3, 6:10)],
                           precursor.retT = c(2:3, 10, 1:5)))
    m <- data.frame(peptide.seq = LETTERS[1:3],
                    precursor.retT.master = 1:3,
                    precursor.retT.slave = c(2:3, 10),
                    deltaRt = c(-1, -1, -7),
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
                precursor.retT = c(0.2, 0.4, 0.6, 1.2, 1.4),
                stringsAsFactors = FALSE)),
             B=list(IdentPeptideData=data.frame(
                peptide.seq = LETTERS[c(1:2, 4)],
                precursor.retT = c(0.21, 0.41, 0.81),
                stringsAsFactors = FALSE)),
             C=list(IdentPeptideData=data.frame(
                peptide.seq = LETTERS[c(1:2, 5)],
                precursor.retT = c(0.22, 0.42, 1.02),
                stringsAsFactors = FALSE)))
   m <- data.frame(peptide.seq = LETTERS[c(1:3, 6:7, 4:5)],
                   precursor.retT = c(0.2, 0.4, 0.6, 1.2, 1.4, 0.8, 1.0),
                   stringsAsFactors = FALSE)
   expect_error(synapter:::.mergeMaster(l[1]),
                "To create a master at least two identification data are needed.")
   expect_equal(suppressWarnings(
                  synapter:::.mergeMaster(l, verbose = FALSE)), m,
                tolerance = 0.02)
})
