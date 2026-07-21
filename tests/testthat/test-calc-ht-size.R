test_that("calc_ht_size returns positive dimensions in requested units", {
  ht <- ComplexHeatmap::Heatmap(
    matrix(seq_len(25), nrow = 5),
    width = grid::unit(30, "mm"),
    height = grid::unit(25, "mm")
  )
  size.in <- calc_ht_size(ht, "inch")
  size.mm <- calc_ht_size(ht, "mm")

  expect_named(size.in, c("width", "height"))
  expect_true(all(is.finite(size.in) & size.in > 0))
  expect_equal(unname(size.mm / size.in), c(25.4, 25.4), tolerance = 1e-5)
})

test_that("calc_ht_size validates inputs", {
  expect_error(calc_ht_size(matrix(1)), "Heatmap or HeatmapList")
  ht <- ComplexHeatmap::Heatmap(matrix(1))
  expect_error(calc_ht_size(ht, NA_character_), "grid unit name")
})
