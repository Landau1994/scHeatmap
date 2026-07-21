# scHeatmap usage gallery

This vignette shows five common single-cell heatmap designs. The examples use
the same interface but make different choices about the heatmap columns and
their ordering.

```r
library(scHeatmap)
library(Seurat)

set.seed(1)
genes <- paste0("gene", 1:12)
cells <- paste0("cell", 1:120)
counts <- matrix(rpois(length(genes) * length(cells), 3),
  nrow = length(genes), dimnames = list(genes, cells))
seu <- CreateSeuratObject(counts)
seu$cell_type <- factor(rep(c("Homeostatic", "Activated", "Cycling"), each = 40))
seu$condition <- rep(rep(c("Control", "Treatment"), each = 20), 3)
seu$sample <- rep(c("sample1", "sample2", "sample3", "sample4"), each = 30)
seu$pseudotime <- ave(seq_len(ncol(seu)), seu$cell_type,
  FUN = function(x) seq(0, 1, length.out = length(x)))
seu <- NormalizeData(seu, verbose = FALSE)

marker_groups <- list(
  Homeostatic = genes[1:4],
  Activated = genes[5:8],
  Cycling = genes[9:12]
)
```

## 1. Standard marker heatmap

Split cells by cell type and preserve their input order.

```r
ht <- sc_heatmap(
  seu,
  features = genes,
  label.features = c("gene1", "gene5", "gene9"),
  split.by = "cell_type",
  feature.groups = marker_groups,
  annotations = c("cell_type", "condition")
)
ComplexHeatmap::draw(ht)
```

All genes remain in the heatmap, while `label.features` marks selected rows
with linked labels. A `cell_type` annotation automatically uses the built-in divergent
categorical palette unless `annotation.colors` supplies an override.

## 2. Balanced downsampled heatmap

Retain the same maximum number of cells from each cell type. `seed` makes the
selection reproducible.

```r
ht <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  downsample = 20,
  seed = 2026,
  annotations = c("cell_type", "condition")
)
```

## 3. Average expression heatmap

Each column represents one cell-type and condition combination instead of one
cell. This design is useful for compact comparisons between biological groups.

```r
ht <- sc_heatmap(
  seu,
  features = genes,
  mode = "average",
  split.by = "cell_type",
  aggregate.by = c("cell_type", "condition"),
  annotations = c("cell_type", "condition"),
  show.column.names = TRUE
)
```

## 4. Continuous or pseudotime ordering

Cells are ordered by a numeric metadata field within each cell-type slice.
`order.by` can also be a feature vector, in which case its mean expression is
used as the ordering score.

```r
ht <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  order.by = "pseudotime",
  annotations = c("condition", "pseudotime")
)

score_ordered <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  order.by = marker_groups$Activated
)
```

## 5. Slice-specific marker ordering

This is the specialized design that motivated scHeatmap. Every experimental
group may use a different marker program and sorting direction.

```r
seu$treatment <- factor(rep(c("WT", "APP", "APP+IL33"), each = 40),
  levels = c("WT", "APP", "APP+IL33"))

ht <- sc_heatmap(
  seu,
  features = genes,
  split.by = "treatment",
  feature.groups = marker_groups,
  sort.by = list(
    WT = marker_groups$Homeostatic,
    APP = marker_groups$Activated,
    `APP+IL33` = marker_groups$Cycling
  ),
  decreasing = c(WT = FALSE, APP = FALSE, `APP+IL33` = TRUE),
  colors = sc_heatmap_palette("purple_black_yellow")
)
```

For an exact project-specific reproduction, see `development/demo.R` in the
source repository.
