context("fragmentlibrary")

test_that(".recalculateFragmentIntensities", {

  fragments <- data.frame(precursor.inten=c(rep(1e3, 3),
                                            rep(2e3, 3),
                                            rep(1e4, 2),
                                            rep(2e4, 2)),
                          product.inten=c(200, 500, 300,
                                          600, 1000, 400,
                                          8e3, 2e3,
                                          15e3, 5e3),
                          fragment.str=c("b1", "b2", "b3",
                                         "b1", "b2", "b4",
                                         "b1", "b4",
                                         "b1", "b2"),
                          peptide.seq=c(rep("ABC", 8),
                                        rep("DEF", 2)),
                          product.rank=c(1:3, 1:3, 1:2, 1:2),
                          run=rep(c(1:3, 1), c(3, 3, 2, 2)),
                          stringsAsFactors=FALSE)

  # round(mean(c(1e3, 2e3, 1e4))) == 4333
  results <- data.frame(precursor.inten=c(rep(4333, 4),
                                          rep(2e4, 2)),
                        product.inten=round(c(8800/13e3*4333,
                                              1500/3e3*4333,
                                              300/1e3*4333,
                                              2e3/1e4*4333,
                                              15e3, 5e3)),
                        fragment.str=c("b1", "b2", "b3", "b4",
                                       "b1", "b2"),
                        peptide.seq=c(rep("ABC", 4),
                                      rep("DEF", 2)),
                        product.rank=c(1:4, 1:2),
                        run=c(rep(1, 3), 2, 1, 1),
                        stringsAsFactors=FALSE)

  expect_equal(synapter:::.recalculateFragmentIntensities(fragments),
               results)
})

