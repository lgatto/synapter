context("utils")

test_that("flatMatchedEMRTs", {
  df <- data.frame(Function=c(1, 1, 2, 2, 1),
                   spectrumID=1:5,
                   precursor.leID.quant=c(1, 3, 3, 9, NA),
                   matched.quant.spectrumIDs=c("1", "3", "3,4", "3,4", ""),
                   stringsAsFactors=FALSE)
  rdf <- data.frame(Function=c(1, 1, 2, 2, 2, 2),
                    spectrumID=c(1:4, 3, 4),
                    precursor.leID.quant=c(1, 3, 3, 3, 9, 9),
                    matched.quant.spectrumIDs=c(1:4, 3, 4),
                    matchType=c("unique-true", "unique-false",
                                "non-unique-true", rep("non-unique-false", 3)),
                    stringsAsFactors=FALSE)

  expect_equal(synapter:::flatMatchedEMRTs(df), rdf)
})

test_that("matched.quant.spectrumIDs2numeric", {
  ids <- c("1", "1,2,3,4", "")
  expect_equal(synapter:::matched.quant.spectrumIDs2numeric(ids),
               list(1, 1:4, numeric(0)))
})
