.sc_heatmap_palettes <- list(
  purple_black_yellow = c("#6B53A0", "#000000", "#F0E926"),
  skyblue_black_orange = c("#87CEEB", "#000000", "#FFA500"),
  coolwarm = c(
    "#3B4CC0", "#6585EC", "#93B5FF", "#C0D4F5", "#E5D8D1",
    "#F7B89C", "#EE8568", "#CC4039", "#B40426"
  ),
  orange_blue = c(
    "#9F3D22", "#DB6525", "#F1AC73", "#D8D4C9", "#A1BCCF",
    "#6FA3CB", "#4171A1", "#2B5B8B"
  ),
  okabe_ito = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#000000"
  ),
  sun = c("#39489F", "#39BBEC", "#F9ED36", "#F38466", "#B81F25"),
  fire = c("#000000", "#380000", "#DF0D00", "#FFB402", "#FFF324", "#FFEE9D"),
  ylgnbu = c("#FFFFD9", "#C7E9B4", "#7FCDBB", "#41B6C4", "#225EA8", "#081D58"),
  monet_sunset = c("#272924", "#323C48", "#67879C", "#7D9390", "#D7695A", "#BA9F84")
)

#' List the built-in scHeatmap palettes
#'
#' @return A character vector of palette names accepted by
#'   [scHeatmap::sc_heatmap_palette()].
#' @export
list_sc_heatmap_palettes <- function() {
  names(.sc_heatmap_palettes)
}

#' Get a built-in scHeatmap color palette
#'
#' Returns curated palettes migrated from the original project utilities. This
#' avoids creating many global palette functions and uses only base R color
#' interpolation.
#'
#' @param name One palette name returned by
#'   [scHeatmap::list_sc_heatmap_palettes()].
#' @param n Number of colors to return. When `NULL`, return the palette's anchor
#'   colors without interpolation.
#' @param reverse Reverse the color order.
#'
#' @return A character vector of hexadecimal colors.
#' @export
#'
#' @examples
#' list_sc_heatmap_palettes()
#' sc_heatmap_palette("purple_black_yellow")
#' sc_heatmap_palette("coolwarm", n = 20)
sc_heatmap_palette <- function(name = "purple_black_yellow", n = NULL,
                               reverse = FALSE) {
  if (!is.character(name) || length(name) != 1L || !name %in% names(.sc_heatmap_palettes)) {
    stop(
      "Unknown palette. Choose one of: ",
      paste(names(.sc_heatmap_palettes), collapse = ", "),
      call. = FALSE
    )
  }
  colors <- .sc_heatmap_palettes[[name]]
  if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1 || n != as.integer(n)) {
      stop("`n` must be NULL or one positive integer.", call. = FALSE)
    }
    colors <- grDevices::colorRampPalette(colors)(as.integer(n))
  }
  if (isTRUE(reverse)) rev(colors) else colors
}
