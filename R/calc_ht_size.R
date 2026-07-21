#' Calculate the rendered size of a ComplexHeatmap
#'
#' `calc_ht_size()` is a convenience wrapper around
#' [ComplexHeatmap::ht_size()]. It converts the returned grid units into a
#' numeric width and height suitable for graphics devices such as
#' [grDevices::pdf()] or [grDevices::png()].
#'
#' For reproducible absolute dimensions, set `width` and `height` (or
#' `heatmap_width` and `heatmap_height`) when constructing the heatmap. A
#' heatmap with fully flexible dimensions may otherwise expand to the active
#' graphics device.
#'
#' @param ht A `Heatmap` or `HeatmapList` object.
#' @param unit Output unit understood by grid, such as `"inch"`, `"mm"`, or
#'   `"cm"`.
#'
#' @return A named numeric vector with elements `width` and `height`.
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(25), nrow = 5)
#' ht <- ComplexHeatmap::Heatmap(
#'   mat,
#'   width = grid::unit(30, "mm"),
#'   height = grid::unit(30, "mm")
#' )
#' calc_ht_size(ht, unit = "inch")
calc_ht_size <- function(ht, unit = "inch") {
  if (!inherits(ht, c("Heatmap", "HeatmapList"))) {
    stop("`ht` must be a ComplexHeatmap Heatmap or HeatmapList object.", call. = FALSE)
  }
  if (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)) {
    stop("`unit` must be one non-empty grid unit name.", call. = FALSE)
  }

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  size <- ComplexHeatmap::ht_size(ht)

  c(
    width = grid::convertWidth(size$width, unit, valueOnly = TRUE),
    height = grid::convertHeight(size$height, unit, valueOnly = TRUE)
  )
}
