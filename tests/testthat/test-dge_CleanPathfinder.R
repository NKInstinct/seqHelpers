pathRes <- system.file("extdata", "pathRes_full.Rds", package = "seqHelpers") |>
  readRDS()

cleanRes <- dge_CleanPathfinder(pathRes$B, geneCount = 2)


# Test Error Handling ----------------------------------------------------------
test_that("Input checking is good", {
  expect_error(dge_GetKeggGeneLists(pathRes))
  expect_error(dge_GetKeggGeneLists(data.frame("A" = runif(10), "B" = runif(10))))
})

# Test Outputs
test_that("Output has expected parameters", {
  expect_true(is(cleanRes, "data.frame"))
  expect_equal(dim(cleanRes), c(3, 9))
  expect_match(cleanRes$ID[[1]], "mmu05203")
})

test_that("Filter Stringency is working", {
  expect_equal(nrow(dge_CleanPathfinder(pathRes$B, geneCount = 100)), 0)
  expect_equal(nrow(dge_CleanPathfinder(pathRes$B, foldChange = 50)), 0)
})
