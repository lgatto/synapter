context("readcsv-functions")

content <- ("Protein Key,Peptide.seq,intensity\nA_B,ABCD,1\nA_D,DCBA,2")

test_that("fetchColnames", {
  expect_equal(synapter:::fetchColnames(textConnection(content)),
               c("Protein Key", "Peptide.seq", "intensity"))
  expect_equal(synapter:::fetchColnames(textConnection(content), sep=" "),
               c("Protein", "Key,Peptide.seq,intensity"))
})

test_that("createReadrColTypes", {
  expect_equal(synapter:::createReadrColTypes(textConnection(content),
                                              keepCols=c("Protein Key" = "c",
                                                         "Peptide.seq" = "c",
                                                         "intensity" = "d")),
               "ccd")
  expect_equal(synapter:::createReadrColTypes(textConnection(content),
                                              keepCols=c("Protein Key" = "c",
                                                         "intensity" = "d")),
               "c_d")
})

test_that("createReadrColTypes", {
  f <- tempfile()
  con <- file(f)
  writeLines(content, con)
  close(con)

  df <- data.frame("ProteinKey"=c("A_B", "A_D"), intensity=1:2,
                   stringsAsFactors=FALSE)
  colnames(df) <- c("Protein Key", "intensity")
  class(df) <- c("tbl_df", "tbl", "data.frame")

  expect_equal(synapter:::readcsv(f, keepCols=c("Protein Key" = "c",
                                                "intensity" = "d")), df)
  unlink(f)
})
