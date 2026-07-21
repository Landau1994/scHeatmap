# scHeatmap

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![简体中文](https://img.shields.io/badge/语言-简体中文-red.svg)](README.zh-CN.md)

`scHeatmap` is an R package for creating publication-ready single-cell
heatmaps directly from Seurat objects. It provides a concise interface to
`ComplexHeatmap` while keeping feature order, marker groups, metadata groups,
and within-group cell ordering explicit and reproducible.

![Example single-cell heatmap](man/figures/CR_2020_GSEreplicate_figS2E_20260721.jpg)

## Features

- Read expression from Seurat v5 assays and layers.
- Draw cell-level or aggregated average-expression heatmaps.
- Split columns using any Seurat metadata field.
- Balance large groups with reproducible per-group downsampling.
- Arrange genes into ordered marker groups.
- Order cells by metadata, pseudotime, custom scores, clustering, or marker sets.
- Add multiple categorical or continuous metadata annotations.
- Apply per-gene z-score scaling, clipping, and custom color schemes.
- Return a standard `ComplexHeatmap` object for further customization.
- Export consistent PDF and PNG figures.

## Requirements

- R >= 4.1
- SeuratObject
- ComplexHeatmap
- circlize

Seurat is recommended for creating and preprocessing input objects.

## Installation

### From GitHub

Once this repository has been pushed to GitHub, install the development version
with:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Landau1994/scHeatmap")
```

Because the example RDS file is stored with Git LFS, install Git LFS before
cloning the full repository:

```bash
git lfs install
git clone https://github.com/Landau1994/scHeatmap.git
```

Installing with `remotes::install_github()` does not require cloning the example
data manually.

### From a local clone

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install("path/to/scHeatmap")
```

## Quick start

```r
library(scHeatmap)

markers <- c("Cd74", "H2-Aa", "H2-Ab1", "Apoe", "Cst7", "Lpl")

ht <- sc_heatmap(
  seu,
  features = markers,
  label.features = c("Cd74", "Apoe", "Lpl"),
  split.by = "cell_type",
  downsample = 100,
  annotations = c("cell_type", "condition")
)

ComplexHeatmap::draw(ht)
save_sc_heatmap(ht, "my_heatmap", width = 6, height = 4.5)
```

For a compact group-level view, aggregate cells by one or more metadata fields:

```r
average_ht <- sc_heatmap(
  seu,
  features = markers,
  mode = "average",
  split.by = "cell_type",
  aggregate.by = c("cell_type", "condition"),
  annotations = c("cell_type", "condition"),
  show.column.names = TRUE
)
```

`order.by` supports metadata such as pseudotime, marker feature sets, custom
numeric scores, and `"cluster"`. The original slice-specific marker ordering is
still available through `sort.by`.

`label.features` keeps every feature in the expression matrix but displays only
selected genes with connecting lines to their heatmap rows. When `cell_type` is included in `annotations`,
its default colors use the original `divergentcolor` palette; supply
`annotation.colors` to override them.

Linked gene labels are italic by default. Use `label.fontface = "plain"` for
upright text, or `label.gp = grid::gpar(fontface = "plain", col = "red")` for
full style control.

Display controls are independent: `show.feature.names` hides gene labels,
`show.column.names` hides cell/aggregate column names, `show.group.names` hides
split titles, and `show.group.annotation` hides the top metadata annotation.
Every heatmap slice has a black border by default. Use `border.color = FALSE`
to remove borders or provide another color to customize them.

By default, expression is read from the active assay's `data` layer, z-scored
per gene, and clipped to `[-2, 2]`. In `mode = "average"`, the arithmetic mean
is first calculated for each gene across cells in every `aggregate.by` group;
the same scaling and clipping are then applied to the group-level matrix. Use
`scale = "none", clip = NULL` to inspect the unscaled group means.

See the [usage gallery](inst/examples/scHeatmap-gallery.md) for five complete
patterns.

See `?sc_heatmap` and `?save_sc_heatmap` for all options.

Calculate device dimensions for a heatmap with fixed body dimensions:

```r
size <- calc_ht_size(ht, unit = "inch")
save_sc_heatmap(ht, "my_heatmap", width = size["width"], height = size["height"])
```

## Built-in palettes

Curated palettes from the original project utilities are available through one
documented interface:

```r
list_sc_heatmap_palettes()

sc_heatmap_palette("purple_black_yellow")
sc_heatmap_palette("coolwarm", n = 20)

ht <- sc_heatmap(
  seu,
  features = markers,
  colors = sc_heatmap_palette("skyblue_black_orange")
)
```

The complete analysis object is stored under `data-raw/` with Git LFS and is
excluded from the installed R package. The original reproducible analysis and
legacy utility definitions are retained under `development/`.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

## License

MIT License.
