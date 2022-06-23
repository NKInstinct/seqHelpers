
data <- system.file("extdata", "dgeRes_full.Rds", package = "seqHelpers") |>
  readRDS()

mat <- dge_GetMatrix(data$DGEList)

gg <- dge_PlotHeatmap(mat)

# Test error handling ----------------------------------------------------------
test_that("non-matrix input is handled", {
  expect_error(dge_PlotHeatmap(data))
})

# Test outputs -----------------------------------------------------------------

test_that("Output is the correct ggplot", {
  expect_true(is(gg, "ggplot"))
})
