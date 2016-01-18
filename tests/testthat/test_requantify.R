context("requantify")

im <- matrix(c(1000, 100, 10, 1, NA,
               800, 100, 10, 1, NA,
               2000, 200, 20, 2, NA,
               6000, 1000, 100, 10, 1), nrow=4, ncol=5, byrow=TRUE,
               dimnames=list(paste("isotopicDistr", paste0("run", 1:4),
                                   sep="."),
                             c("1_0", "1_1", "1_2", "1_3", "1_4")))

test_that(".isUnsaturatedIsotope", {
    munsat <- im
    mode(munsat) <- "logical"
    munsat[] <- TRUE
    munsat[1:3, 5] <- NA
    msat <- munsat
    msat[3:4, 1] <- FALSE
    expect_equal(synapter:::.isUnsaturatedIsotope(im, Inf), munsat)
    expect_equal(synapter:::.isUnsaturatedIsotope(im, 1500), msat)
})

test_that(".runsUnsaturated", {
    iso <- c("1_0", "1_1", "1_2", "1_3", "1_4")
    expect_equal(synapter:::.runsUnsaturated(im, Inf),
                 setNames(rep(TRUE, 5), iso))
    expect_equal(synapter:::.runsUnsaturated(im, 1500),
                 setNames(rep(c(FALSE, TRUE), c(1, 4)), iso))
})

test_that(".isCommonIsotope", {
    expect_equal(synapter:::.isCommonIsotope(im),
                 setNames(rep(c(TRUE, FALSE), c(4, 1)),
                          c("1_0", "1_1", "1_2", "1_3", "1_4")))
})

test_that(".sumIsotopes", {
  m <- matrix(c(1000, 100, 10, NA, NA,
                10000, 1000, NA, 800, 100), byrow=TRUE, nrow=2, ncol=5,
                dimnames=list(c("run1", "run2"),
                              c("1_0", "1_1", "1_2", "2_0", "2_1")))
  r <- matrix(c(1000, 100, 10,
                10800, 1100, 0), byrow=TRUE, nrow=2, ncol=3,
                dimnames=list(c("run1", "run2"),
                              c("0", "1", "2")))
  expect_equal(synapter:::.sumIsotopes(m), r)
})

test_that("requantifySum, all isotopes below threshold", {
  int <- setNames(c(1111, 911, 2222, 7111), colnames(f))
  sat5000_int <- setNames(c(111, 111, 222, 1111), colnames(f))

  expect_equal(synapter:::.requantifySum(im, Inf,
                                         onlyCommonIsotopes=FALSE),
               int)
  expect_equal(synapter:::.requantifySum(im, 5000,
                                         onlyCommonIsotopes=FALSE),
               sat5000_int)
})

test_that("requantifySum, only common isotopes", {
  int <- setNames(c(1111, 911, 2222, 7110), colnames(f))
  sat5000_int <- setNames(c(111, 111, 222, 1110), colnames(f))

  expect_equal(synapter:::.requantifySum(im, Inf,
                                         onlyCommonIsotopes = TRUE),
               int)
  expect_equal(synapter:::.requantifySum(im, 5000,
                                         onlyCommonIsotopes = TRUE),
               sat5000_int)
})

test_that("requantifyReferenceRun", {

  int <- setNames(c(1111, 911, 2222, 7111), colnames(f))
  sat5000_int <- setNames(c(1111, 911, 2222, 11111), colnames(f))

  expect_equal(synapter:::.requantifyReferenceRun(im, Inf), int)
  expect_equal(synapter:::.requantifyReferenceRun(im, 5000), sat5000_int)

  # nothing in common
  notCommon <- matrix(c(1000, 100, 10, NA, NA,
                        NA, NA, NA, 800, 100), byrow=TRUE, nrow=2, ncol=5,
                      dimnames=list(c("run1", "run2"),
                                    c("1_0", "1_1", "1_2", "2_0", "2_1")))
  expect_equal(synapter:::.requantifyReferenceRun(notCommon, Inf),
               setNames(c(1110, 900), c("run1", "run2")))
  expect_equal(synapter:::.requantifyReferenceRun(notCommon, 1000),
               rep(NA_real_, 2))
})