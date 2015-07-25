context("utils")

test_that(".rescueEMRTS", {
  MatchedEMRTs <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                             Counts=1:5,
                             precursor.leID.ident=1:5,
                             idSource="transfer", stringsAsFactors=FALSE)
  MergedEMRTs <- data.frame(precursor.leID.ident=5:1,
                            precursor.inten.quant=10:6)
  rescue <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                       Counts=c(1, 7, 3, 4, 10),
                       precursor.leID.ident=1:5,
                       idSource=c("transfer", "rescue",
                                  "transfer", "transfer", "rescue"),
                       stringsAsFactors=FALSE)
  copy <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                     Counts=6:10,
                     precursor.leID.ident=1:5,
                     idSource="copy", stringsAsFactors=FALSE)
  expect_error(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs, method="foo"))
  expect_equal(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs), rescue)
  expect_equal(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs,
                                       method = "copy"), copy)
})

test_that("flatMatchedEMRTs", {
  df <- data.frame(matchedEMRTs=c(1, 1, 2, 2, 1),
                   spectrumID=1:5,
                   precursor.leID.quant=c(1, 3, 3, 9, NA),
                   matched.quant.spectrumIDs=c("1", "3", "3;4", "3;4", ""),
                   ion_z=c(1, 1, NA, NA, 1),
                   stringsAsFactors=FALSE)
  pep3d <- data.frame(spectrumID=1:4,
                      Function=1:4,
                      ion_z=1:4,
                      stringsAsFactors=FALSE)
  rdf <- data.frame(matchedEMRTs=c(1, 1, 2, 2, 2, 2),
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

test_that(".splitIsotopicDistr", {
  iso <- c("4_0:5;4_1:20;4_2:30",
           "4_0:10;4_1:20;4_4:30;1_1:0",
           "4_0:8;4_2:30")
  r <- list(setNames(c(5, 20, 30), c("4_0", "4_1", "4_2")),
            setNames(c(10, 20, 30, 0), c("4_0", "4_1", "4_4", "1_1")),
            setNames(c(8, 30), c("4_0", "4_2")))
  isoNA <- c("4_0:5;4_1:20;4_2:30", NA)
  rNA <- list(setNames(c(5, 20, 30), c("4_0", "4_1", "4_2")),
              setNames(NA_integer_, NA_character_))

  expect_equal(synapter:::.splitIsotopicDistr(iso), r)
  expect_equal(synapter:::.splitIsotopicDistr(isoNA), rNA)
})
