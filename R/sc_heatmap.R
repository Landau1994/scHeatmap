#' Create a single-cell heatmap from a Seurat object
#'
#' `sc_heatmap()` extracts expression for selected features, optionally scales
#' each feature, splits cells using a metadata column, and returns a
#' [ComplexHeatmap::Heatmap()] object. The input Seurat object is not modified.
#'
#' @param object A Seurat object.
#' @param features Character vector of feature names, in the desired row order.
#' @param group.by Metadata column used to split heatmap columns. Use `NULL` for
#'   no column split. Factor levels determine group order; character metadata is
#'   shown in order of first appearance.
#' @param assay Assay to use. Defaults to [SeuratObject::DefaultAssay()].
#' @param layer Expression layer. Defaults to `"data"`. Use `"counts"` for raw
#'   counts or `"scale.data"` for an existing scaled matrix.
#' @param scale Logical; z-score each feature across the selected cells. This is
#'   done locally and does not change `object`.
#' @param cells Optional character vector of cell names to include.
#' @param feature.groups Optional named list of feature vectors, or a factor
#'   aligned with `features`, used to split rows into marker groups.
#' @param sort.by Optional named list defining within-group cell ordering. Each
#'   name must match a group and each value is a feature vector whose mean is
#'   used as the cell score. A single feature vector applies to every group.
#' @param decreasing Sort high-to-low when `TRUE`. May be one value or a named
#'   logical vector by group.
#' @param colors Three or more colors used for the expression color ramp.
#' @param color.breaks Numeric color breakpoints. Defaults to evenly spaced
#'   values across `clip`.
#' @param clip Length-two numeric range used to clip displayed values. Set to
#'   `NULL` to disable clipping.
#' @param group.colors Optional named colors for a top group annotation.
#' @param show.group.annotation Show the group annotation bar when colors are
#'   supplied.
#' @param cluster.rows,cluster.columns Logical clustering controls.
#' @param show.row.names,show.column.names Logical display controls.
#' @param heatmap.name Legend title.
#' @param use.raster Rasterize the heatmap body. Useful for many cells.
#' @param ... Additional arguments passed to [ComplexHeatmap::Heatmap()].
#'
#' @return A `Heatmap` object. Use [ComplexHeatmap::draw()] to render it or
#'   [scHeatmap::save_sc_heatmap()] to save it.
#' @export
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(60, 4), nrow = 6,
#'     dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10)))
#'   seu <- Seurat::CreateSeuratObject(counts)
#'   seu$condition <- rep(c("control", "treated"), each = 5)
#'   seu <- Seurat::NormalizeData(seu, verbose = FALSE)
#'   ht <- sc_heatmap(seu, paste0("gene", 1:6), group.by = "condition")
#' }
sc_heatmap <- function(object, features, group.by = NULL, assay = NULL,
                       layer = "data", scale = TRUE, cells = NULL,
                       feature.groups = NULL, sort.by = NULL,
                       decreasing = FALSE,
                       colors = c("#6B53A0", "black", "#F0E926"),
                       color.breaks = NULL, clip = c(-2, 2),
                       group.colors = NULL, show.group.annotation = TRUE,
                       cluster.rows = FALSE, cluster.columns = FALSE,
                       show.row.names = TRUE, show.column.names = FALSE,
                       heatmap.name = if (scale) "z-score" else "expression",
                       use.raster = TRUE, ...) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  if (!is.character(features) || !length(features) || anyNA(features)) {
    stop("`features` must be a non-empty character vector without missing values.", call. = FALSE)
  }
  if (anyDuplicated(features)) stop("`features` must not contain duplicates.", call. = FALSE)
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (!assay %in% names(object)) stop("Assay not found: ", assay, call. = FALSE)
  cells <- cells %||% colnames(object)
  unknown.cells <- setdiff(cells, colnames(object))
  if (length(unknown.cells)) stop("Cells not found: ", paste(unknown.cells, collapse = ", "), call. = FALSE)

  mat <- tryCatch(
    SeuratObject::LayerData(object = object, assay = assay, layer = layer),
    error = function(e) stop("Could not read layer `", layer, "` from assay `", assay,
      "`: ", conditionMessage(e), call. = FALSE)
  )
  missing.features <- setdiff(features, rownames(mat))
  if (length(missing.features)) stop("Features not found: ", paste(missing.features, collapse = ", "), call. = FALSE)
  mat <- as.matrix(mat[features, cells, drop = FALSE])
  if (scale) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }

  row.split <- make_feature_groups(feature.groups, features)
  column.split <- make_column_split(object, group.by, cells)
  ord <- order_cells(mat, column.split, sort.by, decreasing)
  mat <- mat[, ord, drop = FALSE]
  if (!is.null(column.split)) column.split <- column.split[ord]
  if (!is.null(clip)) {
    if (!is.numeric(clip) || length(clip) != 2L || clip[1] >= clip[2]) {
      stop("`clip` must be NULL or an increasing numeric vector of length two.", call. = FALSE)
    }
    mat[] <- pmax(clip[1], pmin(clip[2], mat))
  }
  if (is.null(color.breaks)) {
    limits <- clip %||% range(mat, finite = TRUE)
    color.breaks <- seq(limits[1], limits[2], length.out = length(colors))
  }
  if (length(colors) != length(color.breaks) || any(diff(color.breaks) <= 0)) {
    stop("`colors` and increasing `color.breaks` must have equal lengths.", call. = FALSE)
  }
  top.annotation <- make_group_annotation(column.split, group.colors, show.group.annotation)

  ComplexHeatmap::Heatmap(
    mat, name = heatmap.name,
    col = circlize::colorRamp2(color.breaks, colors),
    row_split = row.split, column_split = column.split,
    cluster_rows = cluster.rows, cluster_columns = cluster.columns,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_names = show.row.names, show_column_names = show.column.names,
    top_annotation = top.annotation, use_raster = use.raster, ...
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

make_feature_groups <- function(x, features) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (is.null(names(x)) || any(names(x) == "")) stop("`feature.groups` list must be named.", call. = FALSE)
    mapped <- stats::setNames(
      rep(names(x), lengths(x)),
      unlist(x, use.names = FALSE)
    )
    if (anyDuplicated(names(mapped))) stop("Each feature may occur in only one `feature.groups` group.", call. = FALSE)
    absent <- setdiff(features, names(mapped))
    if (length(absent)) stop("Features missing from `feature.groups`: ", paste(absent, collapse = ", "), call. = FALSE)
    extra <- setdiff(names(mapped), features)
    if (length(extra)) warning("Ignoring features in `feature.groups` not requested in `features`.", call. = FALSE)
    return(factor(unname(mapped[features]), levels = names(x)))
  }
  if (length(x) != length(features)) stop("`feature.groups` must align with `features`.", call. = FALSE)
  if (is.factor(x)) x else factor(x, levels = unique(x))
}

make_column_split <- function(object, group.by, cells) {
  if (is.null(group.by)) return(NULL)
  if (!is.character(group.by) || length(group.by) != 1L) stop("`group.by` must be one metadata column name.", call. = FALSE)
  metadata <- object[[]]
  if (!group.by %in% colnames(metadata)) stop("Metadata column not found: ", group.by, call. = FALSE)
  x <- metadata[[group.by]][match(cells, rownames(metadata))]
  if (anyNA(x)) stop("`group.by` contains missing values for selected cells.", call. = FALSE)
  if (is.factor(x)) droplevels(x) else factor(x, levels = unique(x))
}

order_cells <- function(mat, groups, sort.by, decreasing) {
  if (is.null(sort.by)) return(seq_len(ncol(mat)))
  group.names <- if (is.null(groups)) ".all" else levels(groups)
  if (!is.list(sort.by)) sort.by <- stats::setNames(rep(list(sort.by), length(group.names)), group.names)
  if (is.null(names(sort.by)) || any(names(sort.by) == "")) stop("`sort.by` list must be named by group.", call. = FALSE)
  missing.rules <- setdiff(group.names, names(sort.by))
  if (length(missing.rules)) stop("Missing `sort.by` rules for groups: ", paste(missing.rules, collapse = ", "), call. = FALSE)
  dec <- if (length(decreasing) == 1L) stats::setNames(rep(decreasing, length(group.names)), group.names) else decreasing
  if (is.null(names(dec)) && length(dec) == length(group.names)) names(dec) <- group.names
  if (!all(group.names %in% names(dec))) stop("`decreasing` must supply one value per group.", call. = FALSE)
  unlist(lapply(group.names, function(g) {
    idx <- if (is.null(groups)) seq_len(ncol(mat)) else which(groups == g)
    genes <- sort.by[[g]]
    absent <- setdiff(genes, rownames(mat))
    if (length(absent)) stop("Sort features not present in heatmap: ", paste(absent, collapse = ", "), call. = FALSE)
    score <- colMeans(mat[genes, idx, drop = FALSE])
    idx[order(score, decreasing = isTRUE(dec[[g]]), na.last = TRUE)]
  }), use.names = FALSE)
}

make_group_annotation <- function(groups, colors, show) {
  if (is.null(groups) || is.null(colors) || !show) return(NULL)
  missing.colors <- setdiff(levels(groups), names(colors))
  if (is.null(names(colors)) || length(missing.colors)) {
    stop("`group.colors` must be named and cover every group.", call. = FALSE)
  }
  ComplexHeatmap::HeatmapAnnotation(
    group = groups, col = list(group = colors), show_annotation_name = FALSE
  )
}
