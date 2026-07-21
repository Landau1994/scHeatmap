
library(Seurat)
library(ComplexHeatmap)
library(circlize)

## Run this script from the scHeatmap package root.
devtools::load_all(".")






## Data from GSE147495 processed by by ourself
tmp <- "data-raw/CR_2020_scRNAseq_seu_annoted.rds"
seu <- readRDS(file = tmp)
seu$treatment <- plyr::mapvalues(seu$Group,
  from = c("WT","AD","IL33"),
  to = c("WT","APP","APP+IL33"))

seu$treatment <- factor(seu$treatment,
  levels = c("WT","APP","APP+IL33"))
Idents(seu) <- "cell_type"

## replicate CR 2026 fig S2E

genes <- c(
  ## MHC-II genes
  "Cd74","H2-Aa","H2-Ab1","H2-Eb1",
  ## Homeo 
  "Cap1","Fcrls","Sft2d1","Zfp69",
  ## DAM
  "Nrp1","Lpl","Itgax","Apoe","Cst7")

tmp_select <- c("Homeo","DAM")
seu.sub <- subset(seu,cell_type %in% tmp_select)
table(seu.sub$treatment)



## 1. Row split: map each gene to its marker group (keep manual order)
gene_groups <- factor(
  c(rep("MHC-II", 4), rep("Homeo", 4), rep("DAM", 5)),
  levels = c("MHC-II", "Homeo", "DAM"))
names(gene_groups) <- genes

## 2. Expression matrix (genes x cells), z-scored per gene across cells
DefaultAssay(seu.sub) <- "RNA"
seu.sub <- ScaleData(seu.sub, features = genes)
mat <- as.matrix(GetAssayData(seu.sub, layer = "scale.data"))[genes, ]

## Clip extreme z-scores for a cleaner color range
# tmp.cut <- 1
# mat[mat >  tmp.cut] <-  tmp.cut
# mat[mat < -tmp.cut] <- -tmp.cut

## 3. Column split: treatment per cell (factor keeps WT/APP/APP+IL33 order)
col_split <- seu.sub$treatment[colnames(mat)]

## 4. Color scheme (matches the purple-black-yellow used earlier)
col_fun <- colorRamp2(c(-1, 0, 1),
  c("#6b53a0", "black", "#f0e926"))

## 5. Optional top annotation showing the treatment bar
treat_colors <- c(
  "WT" = "#9b9b9b", 
  "APP" = "#4c8fe2", 
  "APP+IL33" = "#d0051b")
top_anno <- HeatmapAnnotation(
  treatment = col_split,
  col = list(treatment = treat_colors),
  show_annotation_name = FALSE)

## Marker gene sets for the three groups
mhc_genes   <- c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1")
homeo_genes <- c("Cap1", "Fcrls", "Sft2d1", "Zfp69")
dam_genes   <- c("Nrp1", "Lpl", "Itgax", "Apoe", "Cst7")

## Per-slice sorting rule: which score to use, and direction (+1 = low->high, -1 = high->low)
sort_rules <- list(
  "WT"       = list(genes = mhc_genes,   dir =  1),   # WT:       MHC-II score, low  -> high
  "APP"      = list(genes = dam_genes,   dir =  1),   # APP:      DAM score,    low  -> high
  "APP+IL33" = list(genes = homeo_genes, dir = -1))   # APP+IL33: Homeo score,  high -> low

## Build a global column order, sorting each treatment slice independently
ord <- unlist(lapply(levels(col_split), function(tr) {
  cols  <- which(col_split == tr)                          # cell indices in this slice
  rule  <- sort_rules[[tr]]
  score <- colMeans(mat[rule$genes, cols, drop = FALSE])   # per-cell score for this slice
  cols[order(rule$dir * score)]                            # sort within the slice
}))

## Reorder matrix and split factor together (keep them in sync)
mat       <- mat[, ord]
col_split <- col_split[ord]



ph <- Heatmap(mat,
  name            = "z-score",
  col             = col_fun,
  ## rows: fixed marker order, split by group, no clustering
  row_split       = gene_groups,
  cluster_rows    = FALSE,
  cluster_row_slices = FALSE,
  row_title_rot   = 0,
  ## columns: split by treatment, cluster cells within each slice
  column_split    = col_split,
  cluster_columns = F,
  show_column_dend = F,
  cluster_column_slices = FALSE,
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  #top_annotation  = top_anno,
  border          = FALSE,
  use_raster      = TRUE,
  row_names_gp = gpar(fontface = "italic"))   # raster for many cells

#ph
save_sc_heatmap(ph,
  "development/figures/CR_2020_GSEreplicate_figS2E",
  width = 6,
  height = 4.5,
  formats = c("pdf", "png"))
