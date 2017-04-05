context("fragment-functions")

test_that(".filterNeutralLoss", {
  df <- data.frame(foo=1:5,
                   Neutral.LossType=c("None", "None",
                                      "NeutralLoss_NH3",
                                      "NeutralLoss_H2O",
                                      "None"),
                   row.names=as.character(1:5),
                   stringsAsFactors=FALSE)

  dfr <- df[c(1, 2, 5), ]

  expect_equal(synapter:::.filterNeutralLoss(df), dfr)
})

test_that(".filterEmptyFragments", {
  df <- data.frame(foo=1:5,
                   Neutral.LossType=c("", "",
                                      "NeutralLoss_NH3",
                                      "NeutralLoss_H2O",
                                      ""),
                   row.names=as.character(1:5),
                   stringsAsFactors=FALSE)

  dfr <- df[3:4,]

  expect_equal(synapter:::.filterEmptyFragments(df), dfr)
})
