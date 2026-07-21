# scHeatmap

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![简体中文](https://img.shields.io/badge/语言-简体中文-red.svg)](README.zh-CN.md)

`scHeatmap` is an R package for creating publication-ready single-cell
heatmaps directly from Seurat objects. It provides a concise interface to
`ComplexHeatmap` while keeping feature order, marker groups, metadata groups,
and within-group cell ordering explicit and reproducible.

## Features

- Read expression from Seurat v5 assays and layers.
- Split columns using any Seurat metadata field.
- Arrange genes into ordered marker groups.
- Sort cells within each group using marker-set scores.
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

markers <- c(
  "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1",
  "Cap1", "Fcrls", "Sft2d1", "Zfp69",
  "Nrp1", "Lpl", "Itgax", "Apoe", "Cst7"
)

marker_groups <- list(
  `MHC-II` = markers[1:4],
  Homeo = markers[5:8],
  DAM = markers[9:13]
)

ht <- sc_heatmap(
  seu,
  features = markers,
  group.by = "treatment",
  feature.groups = marker_groups,
  group.colors = c(
    WT = "#9b9b9b",
    APP = "#4c8fe2",
    `APP+IL33` = "#d0051b"
  ),
  sort.by = list(
    WT = marker_groups[["MHC-II"]],
    APP = marker_groups$DAM,
    `APP+IL33` = marker_groups$Homeo
  ),
  decreasing = c(WT = FALSE, APP = FALSE, `APP+IL33` = TRUE)
)

ComplexHeatmap::draw(ht)
save_sc_heatmap(
  ht,
  "my_heatmap",
  width = 6,
  height = 4.5,
  formats = c("pdf", "png")
)
```

By default, expression is read from the active assay's `data` layer, z-scored
for each gene across selected cells, and clipped to `[-2, 2]`. Set
`scale = FALSE` and `clip = NULL` to display values from the selected layer
without these transformations.

See `?sc_heatmap` and `?save_sc_heatmap` for all options.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

## License

MIT License.
