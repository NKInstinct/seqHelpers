pathRes_full<- system.file("extdata", "pathRes_full.Rds", package = "seqHelpers") |>
     readRDS()

multiCompLists <- purrr::map(pathRes_full, dge_GetKeggGeneLists)
combinedList <- dge_CollapseKeggLists(multiCompLists)

# Test Error handling ----------------------------------------------------------

test_that("Input checking is working", {
  expect_error(dge_CollapseKeggLists(pathRes_full))
})

# Test outputs -----------------------------------------------------------------
test_that("Output has expected properties", {
  expect_true(is(combinedList[[1]], "character"))
  expect_equal(length(combinedList), 27)
  expect_true("Hepatitis B" %in% names(combinedList))
  expect_match(combinedList[[27]], "Acat1")
})
