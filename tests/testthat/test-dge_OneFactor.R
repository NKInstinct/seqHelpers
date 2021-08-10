# Note - many of these tests are actually testing functions in utils.R, but
# since utils is a dumping ground of small functions, I thought it best to keep
# the tests organized by main function rather than a bunch of mixed tests in
# utils.test.

data <- system.file("extdata", "sampleMatrix.Rds", package = "seqHelpers") |>
  readr::read_rds()
grouplist <- factor(c(rep("A", times = 3),
                      rep("B", times = 3),
                      rep("C", times = 3),
                      rep("D", times = 3)))

# Test that user error handling is working -------------------------------------
test_that("Group mismatch is handled", {
  expect_error(dge_OneFactor(data, Groups = grouplist[c(1:11)]))
})

# Test that the comparison specifications give the right length objects --------
test_that("Default comparison of four groups produces three comps", {
  expect_equal(3, length(dge_OneFactor(data,
                                       Groups = grouplist)))
})

test_that("Specifying two comps returns two comps regardless of factors", {
  expect_equal(2, length(dge_OneFactor(data,
                                       Groups = grouplist,
                                       comps = c("B" = "B-A",
                                                 "D" = "D-C"))))
})

# Test that edgeR is working as intended ---------------------------------------
# Note - will need to change this if the example data ever changes!
res <- dge_OneFactor(data, Groups = grouplist, comps = c("B" = "B-A"))
topGene <- res |>
  dplyr::slice_min(order_by = PValue, n = 1) |>
  dplyr::select(GeneID) |>
  as.character()

test_that("Top gene in comparing B to A is stable", {
  expect_match("ENSMUSG00000000753", topGene)
})
