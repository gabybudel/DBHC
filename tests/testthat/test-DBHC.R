testthat::test_that("errors", {
  testthat::expect_error(
    hmm.clust(sequences = 0),
    "sequences must be a sequence data object of type stslist or must be a data.frame"
  )
})
