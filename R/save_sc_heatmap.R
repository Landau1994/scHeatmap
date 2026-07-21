#' Save a single-cell heatmap
#'
#' Saves the same heatmap to PDF and/or PNG with explicit dimensions. File names
#' are stable by default; set `date = TRUE` to append the current date.
#'
#' @param heatmap A `Heatmap` or `HeatmapList` object.
#' @param filename Output path without a file extension.
#' @param width,height Device dimensions in inches.
#' @param formats Any of `"pdf"` and `"png"`.
#' @param dpi Resolution for PNG output.
#' @param legend.position Heatmap legend position passed to
#'   [ComplexHeatmap::draw()].
#' @param date Append `_YYYYMMDD` to `filename`.
#' @param ... Additional arguments passed to [ComplexHeatmap::draw()].
#'
#' @return Invisibly returns the output file paths.
#' @export
save_sc_heatmap <- function(heatmap, filename, width, height,
                            formats = c("pdf", "png"), dpi = 350,
                            legend.position = "right", date = FALSE, ...) {
  if (!inherits(heatmap, c("Heatmap", "HeatmapList"))) stop("`heatmap` must be a ComplexHeatmap object.", call. = FALSE)
  formats <- match.arg(formats, c("pdf", "png"), several.ok = TRUE)
  stem <- if (date) paste0(filename, "_", format(Sys.Date(), "%Y%m%d")) else filename
  directory <- dirname(stem)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  paths <- paste0(stem, ".", formats)
  for (i in seq_along(formats)) {
    if (formats[i] == "pdf") grDevices::pdf(paths[i], width = width, height = height)
    else grDevices::png(paths[i], width = width, height = height, units = "in", res = dpi)
    tryCatch(
      ComplexHeatmap::draw(heatmap, heatmap_legend_side = legend.position, ...),
      finally = grDevices::dev.off()
    )
  }
  invisible(paths)
}
