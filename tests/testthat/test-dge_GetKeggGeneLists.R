pathRes <- system.file("extdata", "pathRes_full.Rds", package = "seqHelpers") |>
              readRDS()

geneList <- dge_GetKeggGeneLists(pathRes$B)


# Test Error Handling ----------------------------------------------------------
test_that("Input checking is good", {
  expect_error(dge_GetKeggGeneLists(pathRes))
  expect_error(dge_GetKeggGeneLists(data.frame("A" = runif(10), "B" = runif(10))))
})

# Test Output ------------------------------------------------------------------
test_that("Output geneList has expected parameters", {
  expect_equal(length(geneList), 15)
  expect_match(names(geneList)[[1]], "B cell receptor signaling pathway")
  expect_match(geneList[[1]], "Lyn")
})
