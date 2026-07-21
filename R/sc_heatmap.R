#' Create a single-cell heatmap from a Seurat object
#'
#' `sc_heatmap()` creates either a cell-level heatmap or an aggregated average
#' expression heatmap. It supports metadata splits, balanced downsampling,
#' several cell-ordering strategies, multiple annotations, and the specialized
#' per-slice marker ordering used by the original scHeatmap example. The input
#' Seurat object is never modified.
#'
#' @param object A Seurat object.
#' @param features Character vector of feature names, in the desired row order.
#' @param mode Either `"cell"` for one column per cell or `"average"` for one
#'   column per metadata combination.
#' @param group.by Deprecated-compatible alias for `split.by`. Existing code may
#'   continue to use it.
#' @param split.by One metadata column used to split heatmap columns. Factor
#'   levels determine slice order.
#' @param aggregate.by Metadata columns used to calculate average expression
#'   when `mode = "average"`. Defaults to `split.by`.
#' @param assay Assay to use. Defaults to [SeuratObject::DefaultAssay()].
#' @param layer Expression layer. Defaults to `"data"`.
#' @param scale Scaling method: `"row"`, `"column"`, or `"none"`. Logical
#'   `TRUE` and `FALSE` remain supported as aliases for `"row"` and `"none"`.
#' @param cells Optional character vector of cell names to include.
#' @param downsample Maximum number of cells retained per `split.by` group in
#'   cell mode. A single number samples globally when `split.by = NULL`.
#' @param seed Random seed used for reproducible downsampling.
#' @param feature.groups Optional named list of feature vectors, or a factor
#'   aligned with `features`, used to split rows into marker groups.
#' @param label.features Optional subset of `features` whose row names should be
#'   displayed with connecting lines. All features remain in the heatmap.
#'   Defaults to displaying every feature as ordinary row names.
#' @param label.side Side on which linked `label.features` are drawn.
#' @param label.fontface Font face for linked gene labels: `"italic"` (default)
#'   or `"plain"`.
#' @param label.gp Optional complete graphical-parameter override for linked
#'   gene labels. For example, `grid::gpar(fontface = "plain", col = "red")`.
#' @param label.link.gp Graphical parameters for label connecting lines.
#' @param order.by Column ordering rule. Use `"input"`, `"cluster"`, a metadata
#'   column name, a feature vector (ordered by its mean score), or a numeric
#'   vector aligned with selected cells/aggregated columns.
#' @param sort.by Legacy specialized ordering: a named list mapping each
#'   `split.by` group to a different feature vector. It takes precedence over
#'   `order.by` and is supported only in cell mode.
#' @param decreasing Sort high-to-low when `TRUE`. May be one value or a named
#'   logical vector for legacy `sort.by` rules.
#' @param annotations Metadata columns shown as top annotations.
#' @param annotation.colors Optional named list of color mappings, one per
#'   annotation. Discrete mappings should be named character vectors.
#' @param colors Three or more colors used for the expression color ramp.
#' @param color.breaks Numeric color breakpoints. Defaults to evenly spaced
#'   values across `clip`.
#' @param clip Length-two numeric range used to clip displayed values. Set to
#'   `NULL` to disable clipping.
#' @param group.colors Deprecated-compatible color mapping for `split.by`.
#' @param show.group.annotation Show metadata annotations above the heatmap.
#'   Set to `FALSE` to hide annotations even when `annotations` is supplied.
#' @param missing.features How unavailable features are handled: `"error"`,
#'   `"warn"`, or `"drop"`.
#' @param cluster.rows,cluster.columns Logical clustering controls. Setting
#'   `order.by = "cluster"` enables column clustering.
#' @param show.row.names Legacy-compatible control for feature row names.
#' @param show.feature.names Show ordinary feature names or linked
#'   `label.features`. This defaults to `show.row.names` and is the preferred
#'   single-cell-oriented argument.
#' @param show.column.names Show individual cell names or aggregated column
#'   names.
#' @param show.group.names Show row and column split titles created by
#'   `feature.groups` and `split.by`.
#' @param heatmap.name Legend title. Defaults to the displayed expression scale.
#' @param border.color Border color drawn around each heatmap slice. Defaults to
#'   black; set to `FALSE` to disable slice borders. A native ComplexHeatmap
#'   `border` argument supplied through `...` takes precedence.
#' @param use.raster Rasterize the heatmap body. Useful for many cells.
#' @param ... Additional arguments passed to [ComplexHeatmap::Heatmap()].
#'
#' @return A `Heatmap` object. Use [ComplexHeatmap::draw()] to render it or
#'   [scHeatmap::save_sc_heatmap()] to save it.
#' @export
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(120, 4), nrow = 6,
#'     dimnames = list(paste0("gene", 1:6), paste0("cell", 1:20)))
#'   seu <- Seurat::CreateSeuratObject(counts)
#'   seu$condition <- rep(c("control", "treated"), each = 10)
#'   seu <- Seurat::NormalizeData(seu, verbose = FALSE)
#'
#'   cell_ht <- sc_heatmap(seu, paste0("gene", 1:6),
#'     split.by = "condition", downsample = 5)
#'   average_ht <- sc_heatmap(seu, paste0("gene", 1:6),
#'     mode = "average", aggregate.by = "condition")
#' }
sc_heatmap <- function(object, features, mode = c("cell", "average"),
                       group.by = NULL, split.by = NULL, aggregate.by = NULL,
                       assay = NULL, layer = "data", scale = "row",
                       cells = NULL, downsample = NULL, seed = 1,
                       feature.groups = NULL, label.features = NULL,
                       label.side = c("right", "left"),
                       label.fontface = c("italic", "plain"),
                       label.gp = NULL,
                       label.link.gp = grid::gpar(col = "grey50"),
                       order.by = "input",
                       sort.by = NULL, decreasing = FALSE,
                       annotations = NULL, annotation.colors = NULL,
                       colors = c("#6B53A0", "black", "#F0E926"),
                       color.breaks = NULL, clip = c(-2, 2),
                       group.colors = NULL, show.group.annotation = TRUE,
                       missing.features = c("error", "warn", "drop"),
                       cluster.rows = FALSE, cluster.columns = FALSE,
                       show.row.names = TRUE,
                       show.feature.names = show.row.names,
                       show.column.names = FALSE, show.group.names = TRUE,
                       heatmap.name = NULL, border.color = "black",
                       use.raster = TRUE, ...) {
  validate_seurat(object)
  mode <- match.arg(mode)
  missing.features <- match.arg(missing.features)
  label.side <- match.arg(label.side)
  label.fontface <- match.arg(label.fontface)
  label.gp <- label.gp %||% grid::gpar(fontface = label.fontface)
  scale.method <- normalize_scale(scale)
  split.by <- resolve_split_by(group.by, split.by)
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  cells <- validate_cells(object, cells)
  metadata <- object[[]][cells, , drop = FALSE]
  validate_metadata_columns(metadata, unique(c(split.by, aggregate.by, annotations)))

  layer.matrix <- tryCatch(
    SeuratObject::LayerData(object = object, assay = assay, layer = layer),
    error = function(e) stop(
      "Could not read layer `", layer, "` from assay `", assay, "`: ",
      conditionMessage(e), call. = FALSE
    )
  )
  features <- validate_features(features, rownames(layer.matrix), missing.features)
  row.split <- make_feature_groups(feature.groups, features)
  labelled.features <- validate_label_features(features, label.features)

  if (mode == "cell") {
    if (!is.null(aggregate.by)) warning("`aggregate.by` is ignored in cell mode.", call. = FALSE)
    sampled <- downsample_cells(cells, metadata, split.by, downsample, seed)
    cells <- sampled$cells
    metadata <- sampled$metadata
    mat <- layer.matrix[features, cells, drop = FALSE]
    column.metadata <- metadata
    column.split <- metadata_factor(metadata, split.by)
    if (!is.null(sort.by)) {
      ord <- order_cells_legacy(mat, column.split, sort.by, decreasing)
    } else {
      ord <- order_columns(mat, column.metadata, column.split, order.by, decreasing)
    }
  } else {
    if (!is.null(sort.by)) stop("`sort.by` is supported only in cell mode.", call. = FALSE)
    if (!is.null(downsample)) warning("`downsample` is ignored in average mode.", call. = FALSE)
    aggregate.by <- aggregate.by %||% split.by
    if (is.null(aggregate.by) || !length(aggregate.by)) {
      stop("`aggregate.by` is required in average mode when `split.by` is NULL.", call. = FALSE)
    }
    validate_metadata_columns(metadata, aggregate.by)
    if (!is.null(split.by) && !split.by %in% aggregate.by) {
      stop("In average mode, `split.by` must also be included in `aggregate.by`.", call. = FALSE)
    }
    aggregated <- aggregate_expression(
      layer.matrix[features, cells, drop = FALSE], metadata, aggregate.by
    )
    mat <- aggregated$matrix
    column.metadata <- aggregated$metadata
    column.split <- metadata_factor(column.metadata, split.by)
    ord <- order_columns(mat, column.metadata, column.split, order.by, decreasing)
  }

  mat <- as.matrix(mat[, ord, drop = FALSE])
  column.metadata <- column.metadata[ord, , drop = FALSE]
  if (!is.null(column.split)) column.split <- column.split[ord]
  mat <- scale_matrix(mat, scale.method)
  mat <- clip_matrix(mat, clip)
  col.fun <- make_color_function(mat, colors, color.breaks, clip)

  annotations <- resolve_annotations(
    annotations, split.by, group.colors, show.group.annotation
  )
  annotation.colors <- resolve_annotation_colors(
    annotation.colors, split.by, group.colors
  )
  top.annotation <- make_annotations(column.metadata, annotations, annotation.colors)
  if (identical(order.by, "cluster") && is.null(sort.by)) cluster.columns <- TRUE
  heatmap.name <- heatmap.name %||% switch(
    scale.method, row = "row z-score", column = "column z-score", none = "expression"
  )
  label.annotation <- make_label_annotation(
    features, labelled.features, label.side, label.gp, label.link.gp,
    show = show.feature.names
  )
  linked.labels <- !is.null(label.features) && show.feature.names

  heatmap.args <- list(
    mat, name = heatmap.name, col = col.fun,
    row_split = row.split, column_split = column.split,
    cluster_rows = cluster.rows, cluster_columns = cluster.columns,
    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
    show_row_names = show.feature.names && !linked.labels,
    show_column_names = show.column.names,
    left_annotation = if (label.side == "left") label.annotation else NULL,
    right_annotation = if (label.side == "right") label.annotation else NULL,
    top_annotation = top.annotation, use_raster = use.raster
  )
  dots <- list(...)
  if (!"border" %in% names(dots)) dots$border <- border.color
  heatmap.args <- c(heatmap.args, dots)
  if (!show.group.names) {
    heatmap.args["row_title"] <- list(NULL)
    heatmap.args["column_title"] <- list(NULL)
  }
  do.call(ComplexHeatmap::Heatmap, heatmap.args)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

validate_seurat <- function(object) {
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
}

validate_cells <- function(object, cells) {
  cells <- cells %||% colnames(object)
  if (!is.character(cells) || !length(cells) || anyDuplicated(cells)) {
    stop("`cells` must contain unique cell names.", call. = FALSE)
  }
  unknown <- setdiff(cells, colnames(object))
  if (length(unknown)) stop("Cells not found: ", paste(unknown, collapse = ", "), call. = FALSE)
  cells
}

validate_features <- function(features, available, action) {
  if (!is.character(features) || !length(features) || anyNA(features)) {
    stop("`features` must be a non-empty character vector without missing values.", call. = FALSE)
  }
  if (anyDuplicated(features)) stop("`features` must not contain duplicates.", call. = FALSE)
  missing <- setdiff(features, available)
  if (length(missing) && action == "error") {
    stop("Features not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(missing) && action == "warn") {
    warning("Dropping features not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  features <- intersect(features, available)
  if (!length(features)) stop("None of the requested features were found.", call. = FALSE)
  features
}

resolve_split_by <- function(group.by, split.by) {
  if (!is.null(group.by) && !is.null(split.by) && !identical(group.by, split.by)) {
    stop("Use either `group.by` or `split.by`, not different values for both.", call. = FALSE)
  }
  split.by <- split.by %||% group.by
  if (!is.null(split.by) && (!is.character(split.by) || length(split.by) != 1L)) {
    stop("`split.by` must be one metadata column name.", call. = FALSE)
  }
  split.by
}

validate_metadata_columns <- function(metadata, columns) {
  columns <- columns[!is.na(columns)]
  missing <- setdiff(columns, colnames(metadata))
  if (length(missing)) stop("Metadata columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
}

normalize_scale <- function(x) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(if (x) "row" else "none")
  match.arg(x, c("row", "column", "none"))
}

scale_matrix <- function(mat, method) {
  if (method == "row") mat <- t(base::scale(t(mat)))
  if (method == "column") mat <- base::scale(mat)
  mat[!is.finite(mat)] <- 0
  mat
}

clip_matrix <- function(mat, clip) {
  if (is.null(clip)) return(mat)
  if (!is.numeric(clip) || length(clip) != 2L || anyNA(clip) || clip[1] >= clip[2]) {
    stop("`clip` must be NULL or an increasing numeric vector of length two.", call. = FALSE)
  }
  mat[] <- pmax(clip[1], pmin(clip[2], mat))
  mat
}

make_color_function <- function(mat, colors, breaks, clip) {
  if (is.null(breaks)) {
    limits <- clip %||% range(mat, finite = TRUE)
    if (limits[1] == limits[2]) limits <- limits + c(-0.5, 0.5)
    breaks <- seq(limits[1], limits[2], length.out = length(colors))
  }
  if (length(colors) != length(breaks) || any(diff(breaks) <= 0)) {
    stop("`colors` and increasing `color.breaks` must have equal lengths.", call. = FALSE)
  }
  circlize::colorRamp2(breaks, colors)
}

make_feature_groups <- function(x, features) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (is.null(names(x)) || any(names(x) == "")) stop("`feature.groups` list must be named.", call. = FALSE)
    mapped <- stats::setNames(rep(names(x), lengths(x)), unlist(x, use.names = FALSE))
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

validate_label_features <- function(features, label.features) {
  if (is.null(label.features)) return(NULL)
  if (!is.character(label.features) || anyNA(label.features) || anyDuplicated(label.features)) {
    stop("`label.features` must be NULL or a unique character vector.", call. = FALSE)
  }
  unavailable <- setdiff(label.features, features)
  if (length(unavailable)) {
    warning(
      "Ignoring `label.features` not present in `features`: ",
      paste(unavailable, collapse = ", "), call. = FALSE
    )
  }
  intersect(features, label.features)
}

make_label_annotation <- function(features, label.features, side,
                                  label.gp, link.gp, show) {
  if (is.null(label.features) || !length(label.features) || !show) return(NULL)
  at <- match(label.features, features)
  mark <- ComplexHeatmap::anno_mark(
    at = at,
    labels = label.features,
    which = "row",
    side = side,
    labels_gp = label.gp,
    link_gp = link.gp
  )
  ComplexHeatmap::rowAnnotation(
    gene = mark,
    show_annotation_name = FALSE
  )
}

metadata_factor <- function(metadata, column) {
  if (is.null(column)) return(NULL)
  x <- metadata[[column]]
  if (anyNA(x)) stop("`split.by` contains missing values.", call. = FALSE)
  if (is.factor(x)) droplevels(x) else factor(x, levels = unique(x))
}

downsample_cells <- function(cells, metadata, split.by, n, seed) {
  if (is.null(n)) return(list(cells = cells, metadata = metadata))
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n)) {
    stop("`downsample` must be NULL or one positive integer.", call. = FALSE)
  }
  groups <- metadata_factor(metadata, split.by)
  if (is.null(groups)) groups <- factor(rep(".all", length(cells)))
  old.exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old.exists) old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (old.exists) assign(".Random.seed", old.seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  keep <- unlist(lapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (length(idx) > n) sample(idx, n) else idx
  }), use.names = FALSE)
  keep <- sort(keep)
  list(cells = cells[keep], metadata = metadata[keep, , drop = FALSE])
}

aggregate_expression <- function(mat, metadata, aggregate.by) {
  keys <- lapply(metadata[aggregate.by], function(x) {
    if (is.factor(x)) as.character(x) else as.character(x)
  })
  label <- do.call(paste, c(keys, sep = " | "))
  labels <- unique(label)
  result <- vapply(labels, function(x) {
    Matrix::rowMeans(mat[, label == x, drop = FALSE])
  }, numeric(nrow(mat)))
  if (is.null(dim(result))) result <- matrix(result, ncol = 1L)
  rownames(result) <- rownames(mat)
  colnames(result) <- labels
  first <- match(labels, label)
  grouped.metadata <- metadata[first, aggregate.by, drop = FALSE]
  rownames(grouped.metadata) <- labels
  list(matrix = result, metadata = grouped.metadata)
}

order_columns <- function(mat, metadata, groups, order.by, decreasing) {
  if (is.character(order.by) && length(order.by) == 1L && order.by %in% c("input", "cluster")) {
    return(seq_len(ncol(mat)))
  }
  if (is.character(order.by) && length(order.by) == 1L && order.by %in% colnames(metadata)) {
    score <- metadata[[order.by]]
  } else if (is.character(order.by) && all(order.by %in% rownames(mat))) {
    score <- Matrix::colMeans(mat[order.by, , drop = FALSE])
  } else if (is.numeric(order.by) && !is.null(names(order.by)) &&
             all(colnames(mat) %in% names(order.by))) {
    score <- unname(order.by[colnames(mat)])
  } else if (is.numeric(order.by) && length(order.by) == ncol(mat)) {
    score <- order.by
  } else {
    stop(
      "`order.by` must be 'input', 'cluster', a metadata column, available features, ",
      "or a numeric vector aligned with columns.", call. = FALSE
    )
  }
  groups <- groups %||% factor(rep(".all", ncol(mat)))
  unlist(lapply(levels(groups), function(g) {
    idx <- which(groups == g)
    idx[order(score[idx], decreasing = isTRUE(decreasing), na.last = TRUE)]
  }), use.names = FALSE)
}

order_cells_legacy <- function(mat, groups, sort.by, decreasing) {
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
    score <- Matrix::colMeans(mat[genes, idx, drop = FALSE])
    idx[order(score, decreasing = isTRUE(dec[[g]]), na.last = TRUE)]
  }), use.names = FALSE)
}

resolve_annotations <- function(annotations, split.by, group.colors, show) {
  if (!show) return(NULL)
  if (!is.null(group.colors) && show && !is.null(split.by)) annotations <- unique(c(split.by, annotations))
  annotations
}

resolve_annotation_colors <- function(colors, split.by, group.colors) {
  colors <- colors %||% list()
  if (!is.list(colors)) stop("`annotation.colors` must be a named list.", call. = FALSE)
  if (!is.null(group.colors) && !is.null(split.by)) colors[[split.by]] <- group.colors
  colors
}

make_annotations <- function(metadata, annotations, colors) {
  if (is.null(annotations) || !length(annotations)) return(NULL)
  validate_metadata_columns(metadata, annotations)
  annotation.data <- metadata[annotations]
  if ("cell_type" %in% annotations && is.null(colors$cell_type)) {
    values <- annotation.data$cell_type
    levels <- if (is.factor(values)) levels(droplevels(values)) else unique(as.character(values))
    colors$cell_type <- stats::setNames(divergent_colors(length(levels)), levels)
  }
  for (nm in names(colors)) {
    if (!nm %in% annotations) warning("Ignoring colors for unused annotation: ", nm, call. = FALSE)
  }
  ComplexHeatmap::HeatmapAnnotation(
    df = annotation.data,
    col = colors[intersect(names(colors), annotations)],
    show_annotation_name = TRUE
  )
}
