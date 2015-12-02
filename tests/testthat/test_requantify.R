context("requantify")

f <- data.frame(isotopicDistr.run1="1_0:1000;1_1:100;1_2:10;1_3:1;",
                isotopicDistr.run2="1_0:800;1_1:100;1_2:10;1_3:1;",
                isotopicDistr.run3="1_0:2000;1_1:200;1_2:20;1_3:2;",
                isotopicDistr.run4="1_0:6000;1_1:1000;1_2:100;1_3:10;1_4:1",
                stringsAsFactors=FALSE)
im <- .isotopicDistr2matrix(unlist(f))

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

test_that("requantifySum, all isotopes below threshold", {
  int <- setNames(c(1111, 911, 2222, 7111), colnames(f))
  sat5000_int <- setNames(c(111, 111, 222, 1111), colnames(f))

  expect_equal(synapter:::.requantifySum(im, Inf), int)
  expect_equal(synapter:::.requantifySum(im, 5000), sat5000_int)
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
})
