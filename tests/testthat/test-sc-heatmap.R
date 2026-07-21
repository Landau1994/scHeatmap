make_test_object <- function() {
  counts <- Matrix::Matrix(matrix(seq_len(48), nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("c", 1:8))), sparse = TRUE)
  object <- Seurat::CreateSeuratObject(counts)
  object$condition <- factor(rep(c("A", "B"), each = 4), levels = c("B", "A"))
  object$sample <- rep(c("s1", "s2"), times = 4)
  object$pseudotime <- c(4, 1, 3, 2, 8, 5, 7, 6)
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
  expect_error(sc_heatmap(object, "g1", group.by = "absent"), "Metadata columns not found")
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

test_that("cell mode supports balanced reproducible downsampling", {
  object <- make_test_object()
  ht1 <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    downsample = 2, seed = 42, use.raster = FALSE)
  ht2 <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    downsample = 2, seed = 42, use.raster = FALSE)

  expect_equal(ncol(ht1@matrix), 4)
  expect_identical(colnames(ht1@matrix), colnames(ht2@matrix))
  expect_equal(as.integer(table(ht1@matrix_param$column_split)), c(2L, 2L))
})

test_that("metadata and feature-score ordering work within slices", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    order.by = "pseudotime", scale = "none", clip = NULL, use.raster = FALSE)
  expected <- c("c6", "c8", "c7", "c5", "c2", "c4", "c3", "c1")
  expect_identical(colnames(ht@matrix), expected)

  score.ht <- sc_heatmap(object, paste0("g", 1:3), order.by = c("g1", "g2"),
    scale = "none", clip = NULL, use.raster = FALSE)
  expect_s4_class(score.ht, "Heatmap")

  custom.score <- stats::setNames(rev(seq_len(ncol(object))), colnames(object))
  custom.ht <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    downsample = 2, seed = 9, order.by = custom.score, use.raster = FALSE)
  expect_equal(ncol(custom.ht@matrix), 4)
})

test_that("average mode aggregates metadata combinations", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:4), mode = "average",
    split.by = "condition", aggregate.by = c("condition", "sample"),
    annotations = c("condition", "sample"), use.raster = FALSE)

  expect_equal(dim(ht@matrix), c(4, 4))
  expect_setequal(colnames(ht@matrix), c("A | s1", "A | s2", "B | s1", "B | s2"))
  expect_s4_class(ht@top_annotation, "HeatmapAnnotation")
})

test_that("scaling and missing-feature policies are supported", {
  object <- make_test_object()
  expect_warning(
    ht <- sc_heatmap(object, c("g1", "absent", "g2"),
      missing.features = "warn", scale = "column", use.raster = FALSE),
    "Dropping features"
  )
  expect_equal(rownames(ht@matrix), c("g1", "g2"))
  expect_true(all(abs(colMeans(ht@matrix)) < 1e-10))
  expect_silent(sc_heatmap(object, c("g1", "absent"),
    missing.features = "drop", use.raster = FALSE))
})

test_that("multiple annotations and cluster ordering are accepted", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    order.by = "cluster", annotations = c("condition", "sample", "pseudotime"),
    annotation.colors = list(condition = c(A = "red", B = "blue")),
    use.raster = FALSE)
  expect_true(ht@column_dend_param$cluster)
  expect_s4_class(ht@top_annotation, "HeatmapAnnotation")
})

test_that("label.features creates linked labels without dropping rows", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:5),
    label.features = c("g1", "g3", "g5"),
    feature.groups = list(first = c("g1", "g2"), second = c("g3", "g4", "g5")),
    use.raster = FALSE)
  expect_equal(nrow(ht@matrix), 5)
  expect_false(ht@row_names_param$show)
  expect_s4_class(ht@right_annotation, "HeatmapAnnotation")
  expect_named(ht@right_annotation@anno_list, "gene")

  stem <- file.path(tempdir(), "scHeatmap-linked-label-test")
  path <- save_sc_heatmap(ht, stem, 4, 3, formats = "pdf")
  expect_true(file.exists(path))
  expect_gt(file.info(path)$size, 0)

  left.ht <- sc_heatmap(object, paste0("g", 1:3),
    label.features = "g2", label.side = "left", use.raster = FALSE)
  expect_s4_class(left.ht@left_annotation, "HeatmapAnnotation")
  expect_warning(
    sc_heatmap(object, paste0("g", 1:3), label.features = c("g1", "absent"),
      use.raster = FALSE),
    "not present in `features`"
  )
})

test_that("linked gene labels are italic by default and can be plain", {
  object <- make_test_object()
  italic.ht <- sc_heatmap(object, paste0("g", 1:3),
    label.features = "g1", use.raster = FALSE)
  plain.ht <- sc_heatmap(object, paste0("g", 1:3),
    label.features = "g1", label.fontface = "plain", use.raster = FALSE)
  custom.ht <- sc_heatmap(object, paste0("g", 1:3),
    label.features = "g1",
    label.gp = grid::gpar(fontface = "plain", col = "red"),
    use.raster = FALSE)

  italic.gp <- italic.ht@right_annotation@anno_list$gene@fun@var_env$labels_gp
  plain.gp <- plain.ht@right_annotation@anno_list$gene@fun@var_env$labels_gp
  custom.gp <- custom.ht@right_annotation@anno_list$gene@fun@var_env$labels_gp
  expect_identical(unname(italic.gp$font), 3L)
  expect_identical(unname(plain.gp$font), 1L)
  expect_identical(custom.gp$col, "red")
})

test_that("cell_type annotation uses divergent colors by default", {
  object <- make_test_object()
  object$cell_type <- factor(rep(c("T", "B"), each = 4), levels = c("B", "T"))
  ht <- sc_heatmap(object, paste0("g", 1:3), annotations = "cell_type",
    use.raster = FALSE)
  actual <- substr(
    ht@top_annotation@anno_list$cell_type@color_mapping@colors,
    1, 7
  )
  expect_identical(unname(actual), c("#E41A1C", "#377EB8"))
  expect_identical(names(actual), c("B", "T"))

  custom <- c(B = "#111111", T = "#EEEEEE")
  custom.ht <- sc_heatmap(object, paste0("g", 1:3), annotations = "cell_type",
    annotation.colors = list(cell_type = custom), use.raster = FALSE)
  custom.actual <- substr(
    custom.ht@top_annotation@anno_list$cell_type@color_mapping@colors,
    1, 7
  )
  expect_identical(unname(custom.actual), unname(custom))
})

test_that("feature, column, group, and annotation displays can be hidden independently", {
  object <- make_test_object()
  object$cell_type <- factor(rep(c("T", "B"), each = 4))
  ht <- sc_heatmap(
    object,
    paste0("g", 1:4),
    split.by = "cell_type",
    feature.groups = list(first = c("g1", "g2"), second = c("g3", "g4")),
    label.features = c("g1", "g3"),
    annotations = "cell_type",
    show.feature.names = FALSE,
    show.column.names = FALSE,
    show.group.names = FALSE,
    show.group.annotation = FALSE,
    use.raster = FALSE
  )

  expect_false(ht@row_names_param$show)
  expect_false(ht@column_names_param$show)
  expect_null(ht@right_annotation)
  expect_null(ht@top_annotation)
  expect_null(ht@row_title)
  expect_null(ht@column_title)
})

test_that("heatmap slices have configurable black borders by default", {
  object <- make_test_object()
  default.ht <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    use.raster = FALSE)
  no.border.ht <- sc_heatmap(object, paste0("g", 1:3), split.by = "condition",
    border.color = FALSE, use.raster = FALSE)
  native.border.ht <- sc_heatmap(object, paste0("g", 1:3),
    border.color = "blue", border = "red", use.raster = FALSE)

  expect_identical(default.ht@matrix_param$border, "black")
  expect_false(no.border.ht@matrix_param$border)
  expect_identical(native.border.ht@matrix_param$border, "red")
})

test_that("save_sc_heatmap writes requested formats", {
  object <- make_test_object()
  ht <- sc_heatmap(object, paste0("g", 1:3), use.raster = FALSE)
  stem <- file.path(tempdir(), "scHeatmap-test")
  paths <- save_sc_heatmap(ht, stem, 3, 3, formats = "pdf")
  expect_true(file.exists(paths))
})
