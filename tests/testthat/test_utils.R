context("utils")

test_that("flatMatchedEMRTs", {
  df <- data.frame(Function=c(1, 1, 2, 2, 1),
                   spectrumID=1:5,
                   precursor.leID.quant=c(1, 3, 3, 9, NA),
                   matched.quant.spectrumIDs=c("1", "3", "3;4", "3;4", ""),
                   ion_z=c(1, 1, NA, NA, 1),
                   stringsAsFactors=FALSE)
  pep3d <- data.frame(spectrumID=1:4,
                      Function=1:4,
                      ion_z=1:4,
                      stringsAsFactors=FALSE)
  rdf <- data.frame(Function=c(1, 1, 2, 2, 2, 2),
                    spectrumID=c(1:4, 3, 4),
                    precursor.leID.quant=c(1, 3, 3, 3, 9, 9),
                    matched.quant.spectrumIDs=c(1:4, 3, 4),
                    ion_z=c(1, 2, 3, 4, 3, 4),
                    gridSearchResult=c("unique-true", "unique-false",
                                       "non-unique-true",
                                       rep("non-unique-false", 3)),
                    stringsAsFactors=FALSE)

  expect_error(synapter:::flatMatchedEMRTs(df), ".*pep3d.* is missing")
  expect_equal(synapter:::flatMatchedEMRTs(df, pep3d), rdf)
})

test_that("matched.quant.spectrumIDs2numeric", {
  ids <- c("1", "1;2;3;4", "")
  expect_equal(synapter:::matched.quant.spectrumIDs2numeric(ids),
               list(1, 1:4, numeric(0)))
})

