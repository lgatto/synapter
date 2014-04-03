context("spectrum-functions")

test_that(".createMs2SpectrumFromSpectrumXml", {
  m <- c(2, NA, 3, 4, 6, 5)
  i <- c(NA, 11, 12, 13, 15, 14)

  assignments <- data.frame(le_id=1:2,
                            he_id=c("1,2,3,4,6,5", "8,7"),
                            stringsAsFactors=FALSE)
  ms2 <- cbind(HE_ID=1:8, Mass=c(2, NA, 3:8), Intensity=c(NA, 11:17))

  expect_equal(synapter:::.createMs2SpectrumFromSpectrumXml(1, ms2,
                                                            assignments),
               cbind(mass=3:6, intensity=12:15))
})

test_that(".createMs2SpectrumFromFragments", {
  df <- data.frame(precursor.leID=c(rep(1, 6), 2, 2),
                   product.mhp=c(2, NA, 3, 4, 6, 5, 1, 2),
                   product.inten=c(NA, 11, 12, 13, 15, 14, 1, 2),
                   stringsAsFactors=FALSE)
  expect_equal(synapter:::.createMs2SpectrumFromFragments(1, df),
               cbind(mass=3:6, intensity=12:15))
})

test_that(".createSpectrum", {
  m <- c(2, NA, 3, 4, 6, 5)
  i <- c(NA, 11, 12, 13, 15, 14)

  expect_equal(synapter:::.createSpectrum(m, i),
               cbind(mass=3:6, intensity=12:15))
})

