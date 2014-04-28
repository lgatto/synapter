context("fragment-functions")

test_that(".filterNeutralLoss", {
  df <- data.frame(foo=1:5,
                   Neutral.LossType=c("None", "None",
                                      "NeutralLoss_NH3",
                                      "NeutralLoss_H2O",
                                      "None"),
                   row.names=as.character(1:5),
                   stringsAsFactors=FALSE)

  dfr <- data.frame(foo=c(1, 2, 5),
                    Neutral.LossType=rep("None", 3),
                    row.names=as.character(c(1, 2, 5)),
                    stringsAsFactors=FALSE)

  expect_equal(synapter:::.filterNeutralLoss(df), dfr)
})

