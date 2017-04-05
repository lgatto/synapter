context("spectrum-functions")

test_that(".cumsumIntensities", {
  assaydata <- new.env(hash=TRUE, parent=emptyenv(), size=2)
  assign("foo", new("Spectrum2", mz=1:3, intensity=1:3), envir=assaydata)
  assign("bar", new("Spectrum2", mz=1:6, intensity=1:6), envir=assaydata)
  m <- new("MSnExp",
           assayData = assaydata,
           processingData = new("MSnProcess",
                                processing = "Loaded", files="foobar.csv"),
           featureData = new("AnnotatedDataFrame", data=data.frame(spectrum=1:2,
                                                                   row.names=c("bar", "foo"),
                                                                   stringsAsFactors=FALSE)),
           phenoData = new("NAnnotatedDataFrame", data=data.frame(sampleNames=1L,
                                                                  type="spectrum")))
  expect_equal(synapter:::.cumsumIntensities(m), setNames(c(2, 4, 6:9), 1:6))
})

