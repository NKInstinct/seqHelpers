
data <- system.file("extdata", "dgeRes_full.Rds", package = "seqHelpers") |>
  readRDS()

keggRes <- dge_PrepKeggData(data$B)

# Test that error handling is working ------------------------------------------
test_that("Non dataframe input is handled", {
  expect_error(dge_PrepKeggData(data))
})

# Test that outputs are as expected --------------------------------------------

test_that("Ouput is a dataframe with one column and rownames", {
  expect_true(is(keggRes, "data.frame"))
  expect_equal(dim(keggRes), c(87, 1))
  expect_match(rownames(keggRes)[[1]], "Kat2b")
})
