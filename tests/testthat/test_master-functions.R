context("master")

test_that(".orderForMasterModels", {
    l <- list(A=list(IdentPeptideData=data.frame(A=1:2)),
              B=list(IdentPeptideData=data.frame(B=3:10)),
              C=list(IdentPeptideData=data.frame(C=3:5)))
    expect_equal(synapter:::.orderForMasterModels(l),
                 list(c(2, 3, 1), c(3, 1, 2)))
})

