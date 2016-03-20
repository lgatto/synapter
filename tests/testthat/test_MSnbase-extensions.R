context("MSnbase-extensions")

test_that("synapterPlgsAgreement", {
  expect_error(synapterPlgsAgreement(1:3))

  fdata = new("AnnotatedDataFrame",
              data = data.frame(synapterPlgsAgreement.A = c("agree",
                                                            "disagree",
                                                            "no_transfer"),
                                foo.A = 1:3,
                                synapterPlgsAgreement.B = c("agree",
                                                            "agree",
                                                            "disagree"),
                                foo.B = 1:3,
                                foo.C = 1:3,
                                synapterPlgsAgreement.C = c("agree",
                                                            "disagree",
                                                            "disagree"),
                                nIdentified = c(3, 3, 2),
                                nAgree = c(3, 1, 0),
                                nDisagree = c(0, 2, 2),
                                synapterPlgsAgreementRatio = c(1, 1/3, 0/2),
                                row.names = paste0("Pep", 1:3),
                                stringsAsFactors = FALSE))

  eset <- matrix(c(1, 2, NA, 1:3, 1:3), ncol = 3,
                 dimnames = list(paste0("Pep", 1:3), LETTERS[1:3]))
  m <- new("MSnSet",
           exprs = eset,
           processingData = new("MSnProcess",
                                processing = "Coerced from a 'Synapter' object."),
           annotation = "No annotation",
           featureData = fdata[, 1:6])
  r <- new("MSnSet",
           exprs = eset,
           processingData = new("MSnProcess",
                                processing = "Coerced from a 'Synapter' object."),
           annotation = "No annotation",
           featureData = fdata)
  expect_equal(synapterPlgsAgreement(m), r)
})

