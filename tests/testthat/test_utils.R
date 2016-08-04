context("utils")

test_that(".rescueEMRTS", {
  MatchedEMRTs <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                             Counts=1:5,
                             precursor.leID.ident=1:5,
                             isotopicDistr=paste0("iso", 1:5),
                             idSource="transfer", stringsAsFactors=FALSE)
  MergedEMRTs <- data.frame(precursor.leID.ident=5:1,
                            precursor.leID.quant=16:20)
  Pep3D <- data.frame(spectrumID=1:20,
                      Counts=c(rep(0, 15), 10:6),
                      isotopicDistr=paste0("iso", 1:20),
                      stringsAsFactors=FALSE)
  rescue <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                       Counts=c(1, 7, 3, 4, 10),
                       precursor.leID.ident=1:5,
                       isotopicDistr=paste0("iso", c(1, 19, 3:4, 16)),
                       idSource=c("transfer", "rescue",
                                  "transfer", "transfer", "rescue"),
                       stringsAsFactors=FALSE)
  copy <- data.frame(matchedEMRTs=c(1, 2, 1, 1, 2),
                     Counts=6:10,
                     precursor.leID.ident=1:5,
                     isotopicDistr=paste0("iso", c(20:16)),
                     idSource="copy", stringsAsFactors=FALSE)
  expect_error(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs, Pep3D,
                                       method="foo"))
  expect_equal(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs, Pep3D),
               rescue)
  expect_equal(synapter:::.rescueEMRTs(MatchedEMRTs, MergedEMRTs, Pep3D,
               method = "copy"), copy)
})

test_that("flatMatchedEMRTs", {
  df <- data.frame(matchedEMRTs=c(1, 1, 2, 2, 1, 2),
                   spectrumID=1:6,
                   precursor.leID.quant=c(1, 3, 3, 9, NA, NA),
                   matched.quant.spectrumIDs=c("1", "3", "3;4", "3;4", "", "4;5"),
                   ion_z=c(1, 1, NA, NA, 1, 1),
                   stringsAsFactors=FALSE)
  pep3d <- data.frame(spectrumID=1:5,
                      Function=1:5,
                      ion_z=1:5,
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
  rdf2 <- data.frame(matchedEMRTs=c(1, 1, 1, 2, 2, 2, 2, 2, 2),
                     spectrumID=c(1:2, 5, 3:4, 3:4, 4:5),
                     precursor.leID.quant=c(1, 3, NA, 3, 3, 9, 9, NA, NA),
                     matched.quant.spectrumIDs=c(1:2, 5, 3:4, 3:4, 4:5),
                     ion_z=c(1, 2, 5, 3, 4, 3, 4, 4, 5),
                     gridSearchResult=c("unique-true", "unique-false",
                                        "no-quant-id",
                                        "non-unique-true",
                                        rep("non-unique-false", 3),
                                        rep("no-quant-id", 2)),
                     stringsAsFactors=FALSE)

  expect_error(synapter:::flatMatchedEMRTs(df), ".*pep3d.* is missing")
  expect_equal(synapter:::flatMatchedEMRTs(df, pep3d), rdf)
  expect_equal(synapter:::flatMatchedEMRTs(df, pep3d, na.rm=FALSE), rdf2)
})

test_that("matched.quant.spectrumIDs2numeric", {
  ids <- c("1", "1;2;3;4", "")
  expect_equal(synapter:::matched.quant.spectrumIDs2numeric(ids),
               list(1, 1:4, numeric(0)))
})

test_that(".isotopicDistr2matrix", {
  iso <- c(first="4_0:5;4_1:20;4_2:30",
           second="4_0:10;4_1:20;4_4:30;1_1:0",
           third="4_0:8;4_2:30")
  isodf <- data.frame(first="4_0:5;4_1:20;4_2:30",
                      second="4_0:10;4_1:20;4_4:30;1_1:0",
                      third="4_0:8;4_2:30",
                      stringsAsFactors=FALSE)
  isoAllNA <- data.frame(first=NA_character_,
                         second=NA_character_,
                         third=NA_character_,
                         stringsAsFactors=FALSE)

  m <- matrix(c(NA, 5, 20, 30, NA,
                0, 10, 20, NA, 30,
                NA, 8, NA, 30, NA), byrow=TRUE, nrow=3,
              dimnames=list(c("first", "second", "third"),
                            c("1_1", "4_0", "4_1", "4_2", "4_4")))
  isoNA <- c(first="4_0:5;4_1:20;4_2:30", second=NA, third="4_0:2")
  mNA <- matrix(c(5, 20, 30, rep(NA, 3), 2, NA, NA), byrow=TRUE, nrow=3,
                dimnames=list(c("first", "second", "third"),
                              c("4_0", "4_1", "4_2")))
  mAllNA <- matrix(NA_real_, nrow=3,
                   dimnames=list(c("first", "second", "third"), c("1_0")))

  expect_equal(synapter:::.isotopicDistr2matrix(iso), m)
  expect_equal(synapter:::.isotopicDistr2matrix(isodf), m)
  expect_equal(synapter:::.isotopicDistr2matrix(isoNA), mNA)
  expect_equal(synapter:::.isotopicDistr2matrix(isoAllNA), mAllNA)
})

test_that("rawRetTimeModel", {
  retT <- 11:20
  deltaRt <- 1:10*0.2
  span <- 0.75
  expect_equal(loessModel(retT, deltaRt, span)$fitted,
               loess(deltaRt ~ retT, span = span, degree = 1,
                     family = "symmetric", iterations = 4,
                     surface = "direct")$fitted)
})

test_that("modelRetTime", {
  m <- synapter:::modelRetTime(11:20, 1:10*0.2, 0.75)
  expect_identical(names(m), c("lo", "o", "preds", "sd"))
  expect_identical(m$o, 1:10)
  expect_equal(fitted(m$lo), seq(0.2, 2.0, by=0.2))
  expect_equal(unname(m$preds$fit), fitted(m$lo))
  expect_equal(m$preds$df, 6.213695586)
  expect_identical(length(m$sd), 10L)
})

test_that(".commonColnames", {
  x <- data.frame(A=1, B=2, C=3, D=4, E=5)
  y <- data.frame(C=3, D=4, E=5, F=6, G=7, H=8, I=9, J=10)
  z <- data.frame(F=6, G=7, H=8, I=9, J=10)
  expect_identical(synapter:::.commonColnames(x, y), LETTERS[3:5])
  expect_identical(synapter:::.commonColnames(x, y, exclude = LETTERS[3:4]),
                   LETTERS[5])
  expect_identical(synapter:::.commonColnames(x, z), character())
})

test_that(".duplicated2", {
  expect_identical(synapter:::.duplicated2(c(1:4, 1:2)),
                   rep(c(TRUE, FALSE, TRUE), each=2))
})

test_that(".filterNonUniqueIdentMatches", {
  emrts <- data.frame(matchedEMRTs = c(1, 1, 1, 1, 1, 2, 2),
                      spectrumID = c(1, 1, 2, NA, NA, 2, 5),
                      Counts = 1:7)
  r <- emrts
  r$matchedEMRTs[1:2] <- -2
  r$Counts[1:2] <- NA
  r$matchedMultipleIdentEMRTs <- c(rep(TRUE, 2), rep(FALSE, 5))
  expect_equal(synapter:::.filterNonUniqueIdentMatches(emrts), r)
})
