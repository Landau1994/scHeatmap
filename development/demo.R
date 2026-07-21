## scHeatmap development examples
##
## Run this script from the package root. During package development,
## load_all() exposes the current source without reinstalling the package.

library(Seurat)
devtools::load_all(".")


## ---------------------------------------------------------------------------
## 1. Small reproducible examples
## ---------------------------------------------------------------------------

set.seed(1)
genes <- paste0("gene", 1:12)
cells <- paste0("cell", 1:120)
counts <- matrix(
  rpois(length(genes) * length(cells), 3),
  nrow = length(genes),
  dimnames = list(genes, cells)
)
seu <- CreateSeuratObject(counts)
seu$cell_type <- factor(
  rep(c("Homeostatic", "Activated", "Cycling"), each = 40),
  levels = c("Homeostatic", "Activated", "Cycling")
)
seu$condition <- rep(rep(c("Control", "Treatment"), each = 20), 3)
seu$sample <- rep(c("sample1", "sample2", "sample3", "sample4"), each = 30)
seu$pseudotime <- ave(
  seq_len(ncol(seu)), seu$cell_type,
  FUN = function(x) seq(0, 1, length.out = length(x))
)
seu <- NormalizeData(seu, verbose = FALSE)

marker_groups <- list(
  Homeostatic = genes[1:4],
  Activated = genes[5:8],
  Cycling = genes[9:12]
)

## Standard cell-level heatmap. Because the annotation is named `cell_type`,
## its factor levels automatically receive the divergent categorical colors.
ht_standard <- sc_heatmap(
  seu,
  features = genes,
  label.features = c("gene1", "gene5", "gene9"),
  split.by = "cell_type",
  feature.groups = marker_groups,
  annotations = c("cell_type", "condition"),
  colors = sc_heatmap_palette("coolwarm")
)
ComplexHeatmap::draw(ht_standard, merge_legend = TRUE)

## Explicit colors override the automatic cell_type palette.
cell_type_colors <- c(
  Homeostatic = "#4DAF4A",
  Activated = "#E41A1C",
  Cycling = "#377EB8"
)
ht_custom_annotation <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  annotations = c("cell_type", "condition"),
  annotation.colors = list(cell_type = cell_type_colors),
  colors = sc_heatmap_palette("rdbu")
)

## Balanced downsampling retains at most 20 cells per cell type.
ht_downsampled <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  downsample = 20,
  seed = 2026,
  annotations = c("cell_type", "condition"),
  colors = sc_heatmap_palette("rdbu")
)

## Average expression: one column per cell_type and condition combination.
ht_average <- sc_heatmap(
  seu,
  features = genes,
  mode = "average",
  split.by = "cell_type",
  aggregate.by = c("cell_type", "condition"),
  annotations = c("cell_type", "condition"),
  show.column.names = TRUE
)

## Continuous metadata and marker-score ordering.
ht_pseudotime <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  order.by = "pseudotime",
  annotations = c("condition", "pseudotime")
)
ht_marker_score <- sc_heatmap(
  seu,
  features = genes,
  split.by = "cell_type",
  order.by = marker_groups$Activated
)

## Slice-specific ordering: each treatment uses a different marker program.
seu$treatment <- factor(
  rep(c("WT", "APP", "APP+IL33"), each = 40),
  levels = c("WT", "APP", "APP+IL33")
)
ht_slice_specific <- sc_heatmap(
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


## ---------------------------------------------------------------------------
## 2. PBMC3K marker heatmaps
## ---------------------------------------------------------------------------

if (!requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
  stop("Install pbmc3k.SeuratData before running the PBMC3K examples.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install dplyr before running the PBMC3K examples.")
}

data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- Seurat::UpdateSeuratObject(pbmc3k.final)
pbmc$cell_type <- factor(Seurat::Idents(pbmc), levels = levels(Seurat::Idents(pbmc)))

pbmc_markers <- FindAllMarkers(pbmc, only.pos = TRUE, verbose = FALSE)

## Keep the strongest markers and assign a duplicated marker to only its first
## (highest-ranked) cluster, so feature.groups remains unambiguous.
markers_per_cell_type <- 20
top_marker_table <- pbmc_markers |>
  dplyr::filter(avg_log2FC > 1) |>
  dplyr::group_by(cluster) |>
  dplyr::slice_head(
    n = markers_per_cell_type) |>
  dplyr::ungroup() |>
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) |>
  dplyr::distinct(gene, .keep_all = TRUE)

pbmc_feature_groups <- split(top_marker_table$gene, top_marker_table$cluster)
pbmc_features <- unlist(pbmc_feature_groups, use.names = FALSE)

## Only these canonical markers receive visible row labels. Every selected top
## marker remains in the matrix and contributes to the heatmap.
pbmc_label_genes <- c(
  "IL7R", "CCR7", "CD14", "LYZ", "S100A4", "MS4A1", "CD8A",
  "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "GP9","PF4"
)
pbmc_labels <- intersect(pbmc_label_genes, pbmc_features)

## Full PBMC3K cell-level heatmap with default divergent cell-type colors.
pbmc_ht <- sc_heatmap(
  pbmc,
  features = pbmc_features,
  label.features = pbmc_labels,
  label.fontface = "italic",
  split.by = "cell_type",
  feature.groups = pbmc_feature_groups,
  annotations = "cell_type",
  colors = sc_heatmap_palette("solarExtra"),
  show.feature.names = TRUE,
  show.group.names =  FALSE,
  show.group.annotation = TRUE,
  show.column.names = FALSE,
  border.color = "black"
)
#ComplexHeatmap::draw(pbmc_ht, merge_legend = TRUE)
save_sc_heatmap(
  pbmc_ht,
  "development/figures/pbmc_ht",
  width = 6,
  height = 4.5,
  formats = c("pdf", "png"),
  merge_legend = TRUE
)
## A more compact and balanced version for routine plotting.
pbmc_downsampled_ht <- sc_heatmap(
  pbmc,
  features = pbmc_features,
  label.features = pbmc_labels,
  split.by = "cell_type",
  feature.groups = pbmc_feature_groups,
  downsample = 30,
  seed = 2026,
  annotations = "cell_type",
  colors = sc_heatmap_palette("rdbu"),
  show.feature.names = TRUE,
  show.group.names =  FALSE,
  show.group.annotation = TRUE,
  show.column.names = FALSE,
  border.color = "black"
)
ComplexHeatmap::draw(pbmc_downsampled_ht, merge_legend = TRUE)

## One average-expression column per annotated PBMC cell type.
pbmc_average_ht <- sc_heatmap(
  pbmc,
  features = pbmc_features,
  label.features = pbmc_labels,
  mode = "average",
  split.by = "cell_type",
  aggregate.by = "cell_type",
  feature.groups = pbmc_feature_groups,
  annotations = "cell_type",
  show.column.names = TRUE,
  colors = sc_heatmap_palette("rdbu"),
  show.feature.names = TRUE,
  show.group.names =  FALSE,
  show.group.annotation = TRUE
)
ComplexHeatmap::draw(pbmc_average_ht, merge_legend = TRUE)

## Uncomment to save the balanced PBMC3K heatmap.
# save_sc_heatmap(
#   pbmc_downsampled_ht,
#   "development/figures/pbmc3k_marker_heatmap",
#   width = 9,
#   height = 8,
#   formats = c("pdf", "png")
# )


## ---------------------------------------------------------------------------
## 3. CR2020 Figure S2E-style reproduction
## ---------------------------------------------------------------------------

cr_seu <- readRDS("data-raw/CR_2020_scRNAseq_seu_annoted.rds")
cr_seu$treatment <- factor(
  plyr::mapvalues(
    cr_seu$Group,
    from = c("WT", "AD", "IL33"),
    to = c("WT", "APP", "APP+IL33")
  ),
  levels = c("WT", "APP", "APP+IL33")
)
Seurat::Idents(cr_seu) <- "cell_type"
cr_subset <- subset(cr_seu, cell_type %in% c("Homeo", "DAM"))
SeuratObject::DefaultAssay(cr_subset) <- "RNA"

cr_marker_groups <- list(
  `MHC-II` = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"),
  Homeo = c("Cap1", "Fcrls", "Sft2d1", "Zfp69"),
  DAM = c("Nrp1", "Lpl", "Itgax", "Apoe", "Cst7")
)
cr_features <- unlist(cr_marker_groups, use.names = FALSE)

cr_ht <- sc_heatmap(
  cr_subset,
  features = cr_features,
  split.by = "treatment",
  feature.groups = cr_marker_groups,
  sort.by = list(
    WT = cr_marker_groups[["MHC-II"]],
    APP = cr_marker_groups$DAM,
    `APP+IL33` = cr_marker_groups$Homeo
  ),
  decreasing = c(WT = FALSE, APP = FALSE, `APP+IL33` = TRUE),
  colors = sc_heatmap_palette("purple_black_yellow"),
  color.breaks = c(-1, 0, 1),
  clip = c(-1, 1),
  row_title_rot = 0,
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  border = FALSE,
  row_names_gp = grid::gpar(fontface = "italic")
)

save_sc_heatmap(
  cr_ht,
  "development/figures/CR_2020_GSEreplicate_figS2E",
  width = 6,
  height = 4.5,
  formats = c("pdf", "png")
)
