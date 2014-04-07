context("MSnbase-extensions")

test_that("normalize,Spectrum2", {
  s1 <- new("Spectrum2", mz=1:5, intensity=1:5, precursorIntensity=10)

  ## max is default
  expect_equal(intensity(normalize(s1)), (1:5)/5)
  expect_equal(intensity(normalize(s1, method="sum")), (1:5)/15)
  expect_equal(intensity(normalize(s1, method="precursor")), (1:5)/10)
})

test_that(".commonPeaks", {
  m1 <- new("Spectrum2", mz=1:5, intensity=1:5)
  m2 <- new("Spectrum2", mz=c(1, 2.1, 2.9, 4.5, 4.8), intensity=1:5)
  m3 <- new("Spectrum2", mz=3, intensity=1)

  expect_equal(synapter:::.commonPeaks(m1, m1), rep(TRUE, 5))
  expect_equal(synapter:::.commonPeaks(m1, m2), c(TRUE, rep(FALSE, 4)))
  expect_equal(synapter:::.commonPeaks(m1, m2, tolerance=0.05),
               c(TRUE, FALSE, TRUE, FALSE, TRUE))
  expect_equal(synapter:::.commonPeaks(m1, m3), TRUE)
  expect_equal(synapter:::.commonPeaks(m3, m1),
               c(FALSE, FALSE, TRUE, FALSE, FALSE))
})

test_that(".nCommonPeaks", {
  m1 <- new("Spectrum2", mz=1:5, intensity=1:5)
  m2 <- new("Spectrum2", mz=c(1, 2.1, 2.9, 4.5, 4.8), intensity=1:5)
  m3 <- new("Spectrum2", mz=3, intensity=1)

  expect_true(synapter:::.nCommonPeaks(m1, m1) == 5)
  expect_true(synapter:::.nCommonPeaks(m2, m2) == 5)
  expect_true(synapter:::.nCommonPeaks(m1, m2) == 1)
  expect_true(synapter:::.nCommonPeaks(m1, m2, tolerance=0.04) == 2)
  expect_true(synapter:::.nCommonPeaks(m1, m2, tolerance=0.05) == 3)
  expect_true(synapter:::.nCommonPeaks(m1, m2, tolerance=0.06) == 4)
  expect_true(synapter:::.nCommonPeaks(m1, m2, tolerance=0.2) == 5)
  expect_true(synapter:::.nCommonPeaks(m3, m1) == 1)
})


