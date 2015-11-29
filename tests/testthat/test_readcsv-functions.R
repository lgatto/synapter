context("readcsv-functions")

test_that("readcsv", {
  content <- "Protein Key,Peptide.seq,intensity\nA_B,ABCD,1\nA_D,DCBA,2"

  df <- data.frame("ProteinKey"=c("A_B", "A_D"), intensity=1:2,
                   stringsAsFactors=FALSE)
  colnames(df) <- c("Protein Key", "intensity")

  expect_equal(synapter:::readcsv(content, keepCols=c("Protein Key" = "c",
                                                      "intensity" = "d")), df)
})
