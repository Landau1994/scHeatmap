make_test_object <- function() {
  counts <- Matrix::Matrix(matrix(seq_len(48), nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("c", 1:8))), sparse = TRUE)
  object <- Seurat::CreateSeuratObject(counts)
  object$condition <- factor(rep(c("A", "B"), each = 4), levels = c("B", "A"))
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("sc_heatmap creates a heatmap without changing its input", {
  object <- make_test_object()
  original <- object
  ht <- sc_heatmap(object, paste0("g", 1:4), group.by = "condition",
    feature.groups = list(first = c("g1", "g2"), second = c("g3", "g4")),
    group.colors = c(A = "red", B = "blue"), use.raster = FALSE)
  expect_s4_class(ht, "Heatmap")
  expect_identical(object, original)
})

test_that("sc_heatmap validates features and metadata", {
  object <- make_test_object()
  expect_error(sc_heatmap(object, "absent"), "Features not found")
  expect_error(sc_heatmap(object, "g1", group.by = "absent"), "Metadata column not found")
  expect_error(sc_heatmap(object, c("g1", "g2"),
    feature.groups = list(one = "g1")), "missing from `feature.groups`")
})

test_that("within-group sorting rules are accepted", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:3), group.by = "condition",
    sort.by = list(B = c("g1", "g2"), A = "g3"),
    decreasing = c(B = FALSE, A = TRUE), use.raster = FALSE)
  expect_s4_class(ht, "Heatmap")
})

test_that("save_sc_heatmap writes requested formats", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:3), use.raster = FALSE)
  stem <- file.path(tempdir(), "scHeatmap-test")
  paths <- save_sc_heatmap(ht, stem, 3, 3, formats = "pdf")
  expect_true(file.exists(paths))
})
