context("readcsv-functions")

content <- "Protein Key,Peptide.seq,intensity\nA_B,ABCD,1\nA_D,DCBA,2"
contentQuoted <- "\"Protein Key\",\"Peptide.seq\",\"intensity\"\nA_B,ABCD,1\nA_D,DCBA,2"

test_that("fetchColnames", {
  expect_equal(synapter:::fetchColnames(textConnection(content)),
               c("Protein Key", "Peptide.seq", "intensity"))
  expect_equal(synapter:::fetchColnames(textConnection(content), sep=" "),
               c("Protein", "Key,Peptide.seq,intensity"))
  ## quoted headers; see #96
  expect_equal(synapter:::fetchColnames(textConnection(contentQuoted)),
               c("Protein Key", "Peptide.seq", "intensity"))
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
  ## quoted headers; see #96
  expect_equal(synapter:::createReadrColTypes(textConnection(contentQuoted),
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

  expect_equal(synapter:::readcsv(f, keepCols=c("Protein Key" = "c",
                                                "intensity" = "d")), df)
  unlink(f)
})