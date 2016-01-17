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

