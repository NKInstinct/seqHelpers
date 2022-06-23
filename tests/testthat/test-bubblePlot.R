pathRes <- system.file("extdata", "pathRes_full.Rds", package = "seqHelpers") |>
  readRDS()

plot <- bubblePlot(pathRes$B)


# Test Error Handling ----------------------------------------------------------
test_that("Input checking is good", {
  expect_error(bubblePlot(pathRes))
  expect_error(bubblePlot(data.frame("A" = runif(10), "B" = runif(10))))
})

# Test output ------------------------------------------------------------------
test_that("Output is a ggplot object", {
  expect_true(is(plot, "ggplot"))
})
