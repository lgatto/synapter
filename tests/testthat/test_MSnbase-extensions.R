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

test_that(".correctIntensity", {
  fdata = new("AnnotatedDataFrame",
              data = data.frame(intensityCorrectionFactor.A = 1:3,
                                intensityCorrectionFactor.B = 2:4,
                                intensityCorrectionFactor.C = 3:5,
                                row.names = paste0("Pep", 1:3),
                                stringsAsFactors = FALSE))

  eset <- matrix(c(1, 2, NA, 1:3, 1:3), ncol = 3,
                 dimnames = list(paste0("Pep", 1:3), LETTERS[1:3]))
  m <- new("MSnSet",
           exprs = eset,
           processingData = new("MSnProcess",
                                processing = "Coerced from a 'Synapter' object."),
           annotation = "No annotation",
           featureData = fdata)
  r <- new("MSnSet",
           exprs = eset * c(1:3, 2:4, 3:5),
           processingData = new("MSnProcess",
                                processing = "Coerced from a 'Synapter' object."),
           annotation = "No annotation",
           featureData = fdata)
  expect_equal(synapter:::.correctIntensity(m, method="correct"), r)
  expect_equal(synapter:::.correctIntensity(r, method="undo"), m)
})

test_that(".getSpectrum", {
  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=2)
  assign("foo", new("Spectrum2", mz=1:3, intensity=1:3), envir=assaydata)
  assign("bar", new("Spectrum2", mz=4:6, intensity=4:6), envir=assaydata)
  m <- new("MSnExp",
           assayData = assaydata,
           processingData = new("MSnProcess",
                                processing = "Loaded", files="foobar.csv"),
           featureData = new("AnnotatedDataFrame", data=data.frame(spectrum=1:2,
                                                                   row.names=c("bar", "foo"),
                                                                   stringsAsFactors=FALSE)),
           phenoData = new("NAnnotatedDataFrame", data=data.frame(sampleNames=1L,
                                                                  type="spectrum")))
  expect_equal(synapter:::.getSpectrum(key="spectrum:200", m),
               new("Spectrum2", precScanNum=200L, fromFile=1L,
                   acquisitionNum=NA_integer_,
                   collisionEnergy=NA_real_,
                   precursorCharge=NA_integer_,
                   precursorIntensity=NA_real_,
                   precursorMz=NA_real_,
                   peaksCount=integer(),
                   scanIndex=integer()))
  expect_equal(synapter:::.getSpectrum(key="foo", m), m[["foo"]])
  expect_equal(synapter:::.getSpectrum(key="bar", m), m[["bar"]])
})

test_that(".sumAllPeakCounts", {
  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=2)
  assign("foo", new("Spectrum2", mz=1:3, intensity=1:3), envir=assaydata)
  assign("bar", new("Spectrum2", mz=4:6, intensity=4:6), envir=assaydata)
  m <- new("MSnExp",
           assayData = assaydata,
           processingData = new("MSnProcess",
                                processing = "Loaded", files="foobar.csv"),
           featureData = new("AnnotatedDataFrame", data=data.frame(spectrum=1:2,
                                                                   row.names=c("bar", "foo"),
                                                                   stringsAsFactors=FALSE)),
           phenoData = new("NAnnotatedDataFrame", data=data.frame(sampleNames=1L,
                                                                  type="spectrum")))
  expect_equal(synapter:::.sumAllPeakCounts(m), 6)
})

