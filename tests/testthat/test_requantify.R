context("requantify")

test_that("requantifyReferenceRun", {
  f <- data.frame(isotopicDistr.run1="1_0:1000;1_1:100;1_2:10;1_3:1;",
                  isotopicDistr.run2="1_0:800;1_1:100;1_2:10;1_3:1;",
                  isotopicDistr.run3="1_0:2000;1_1:200;1_2:20;1_3:2;",
                  isotopicDistr.run4="1_0:6000;1_1:1000;1_2:100;1_3:10;1_4:1",
                  stringsAsFactors=FALSE)

  int <- setNames(c(1111, 911, 2222, 7111), colnames(f))
  sat5000_int <- setNames(c(1111, 911, 2222, 11110), colnames(f))

  expect_equal(synapter:::.requantifyReferenceRun(f, Inf), int)
  expect_equal(synapter:::.requantifyReferenceRun(f, 5000), sat5000_int)
})
