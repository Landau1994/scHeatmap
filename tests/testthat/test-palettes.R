test_that("built-in palettes are valid and discoverable", {
  palettes <- list_sc_heatmap_palettes()
  expect_true(all(c("purple_black_yellow", "coolwarm", "okabe_ito", "divergent") %in% palettes))
  expect_true(all(grDevices::col2rgb(sc_heatmap_palette("coolwarm", 12)) >= 0))
  expect_length(sc_heatmap_palette("coolwarm", 12), 12)
})

test_that("legacy divergent colors preserve their original order", {
  expect_identical(
    scHeatmap:::divergent_colors(3),
    c("#E41A1C", "#377EB8", "#4DAF4A")
  )
  expect_length(scHeatmap:::divergent_colors(30), 30)
})

test_that("palette direction and input validation work", {
  colors <- sc_heatmap_palette("sun")
  expect_identical(sc_heatmap_palette("sun", reverse = TRUE), rev(colors))
  expect_error(sc_heatmap_palette("missing"), "Unknown palette")
  expect_error(sc_heatmap_palette("sun", 0), "positive integer")
})
