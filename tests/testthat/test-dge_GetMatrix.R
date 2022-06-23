# Note - many of these tests are actually testing functions in utils.R, but
# since utils is a dumping ground of small functions, I thought it best to keep
# the tests organized by main function rather than a bunch of mixed tests in
# utils.test.

data <- system.file("extdata", "dgeRes_full.Rds", package = "seqHelpers") |>
  readRDS()

mat <- dge_GetMatrix(data$DGEList)

scaleCheck <- list(
  "none" = mat[1,1],
  "rows" = dge_GetMatrix(data$DGEList, scale = "rows")[1,1],
  "cols" = dge_GetMatrix(data$DGEList, scale = "cols")[1,1]
)


# Test that user error handling is working -------------------------------------
test_that("Class checks are in place", {
  expect_error(dge_GetMatrix(c(1, 3, 5, 7, 9)))
})

# Test that the output matrix is correct ---------------------------------------
test_that("Returned object is a matrix", {
  expect_true(is.matrix(mat))
})

test_that("Matrix dimensions are correct", {
  expect_equal(dim(mat), c(308, 12))
})

test_that("First gene name is present & stable", {
  expect_match(rownames(mat)[[1]], "ENSMUSG00000000295")
})

test_that("Sample values are distinct", {
  expect_false(mat[1,1] == mat[1,2])
})

# Test that optional args work -------------------------------------------------

test_that("Scaling works as expected", {
  expect_equal(round(scaleCheck$none, digits = 0), 158)
  expect_equal(round(scaleCheck$rows, digits = 2), -1.14)
  expect_equal(round(scaleCheck$cols, digits = 3), -0.197)
})
