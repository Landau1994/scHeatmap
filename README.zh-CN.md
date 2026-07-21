# scHeatmap

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![简体中文](https://img.shields.io/badge/语言-简体中文-red.svg)](README.zh-CN.md)

`scHeatmap` 是一个从 Seurat 对象直接创建单细胞热图的 R 包。它在
`ComplexHeatmap` 基础上提供简洁接口，同时明确保留基因顺序、marker 分组、
元数据分组以及组内细胞排序，使绘图过程更容易复现和调整。

![使用 scHeatmap 复现的文献风格小胶质细胞热图](man/figures/CR_2020_GSEreplicate_figS2E_20260721.jpg)

*该文献风格小胶质细胞热图使用 Lau 等（2020）的研究数据复现，数据来源论文为
[《IL-33-PU.1 Transcriptome Reprogramming Drives Functional State Transition
and Clearance Activity of Microglia in Alzheimer’s Disease》](https://www.cell.com/cell-reports/fulltext/S2211-1247(20)30430-7)。
完整代码见[文献图复现](#文献图复现)。*

## 主要功能

- 支持读取 Seurat v5 assay 和 layer。
- 支持逐细胞热图和聚合平均表达热图。
- 使用任意 Seurat metadata 字段对细胞分组。
- 对不同分组进行可复现的等量下采样。
- 按指定顺序组织基因和 marker 集合。
- 按 metadata、pseudotime、自定义分数、聚类或 marker 集合排序。
- 添加多个分类或连续 metadata 注释。
- 支持逐基因 z-score、数值截断和自定义配色。
- 返回标准 `ComplexHeatmap` 对象，方便继续添加注释或调整样式。
- 将同一热图稳定导出为 PDF 和 PNG。

## 环境要求

- R >= 4.1
- SeuratObject
- ComplexHeatmap
- circlize

推荐使用 Seurat 创建和预处理输入对象。

## 安装

### 从 GitHub 安装

仓库上传到 GitHub 后，可以安装开发版本：

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Landau1994/scHeatmap")
```

仓库中的示例 RDS 文件使用 Git LFS 管理；如果需要克隆完整仓库，应先安装
Git LFS：

```bash
git lfs install
git clone https://github.com/Landau1994/scHeatmap.git
```

直接使用 `remotes::install_github()` 安装 R 包时，不需要手工克隆示例数据。

### 从本地仓库安装

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install("path/to/scHeatmap")
```

## 快速开始

下面使用 SeuratData 提供的 PBMC3K 对象给出可复现示例。如果尚未安装该数据集，
先运行一次 `SeuratData::InstallData("pbmc3k")`。使用 `FindAllMarkers()` 自动筛选
marker 的完整流程保留在 [`development/demo.R`](development/demo.R) 中。

```r
library(scHeatmap)

data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- Seurat::UpdateSeuratObject(pbmc3k.final)
pbmc$cell_type <- factor(
  Seurat::Idents(pbmc),
  levels = levels(Seurat::Idents(pbmc))
)

marker_groups <- list(
  `T cells` = c("IL7R", "CCR7", "LTB", "CD3D", "CD8A", "S100A4"),
  `B cells` = c("MS4A1", "CD79A", "CD37", "CD79B"),
  Monocytes = c("CD14", "LYZ", "FCGR3A", "MS4A7", "CTSS"),
  NK = c("GNLY", "NKG7"),
  Dendritic = c("FCER1A", "CST3"),
  Platelets = c("GP9", "PF4")
)
markers <- unlist(marker_groups, use.names = FALSE)
label_genes <- c("IL7R", "CD8A", "MS4A1", "CD14", "FCGR3A", "GNLY", "FCER1A", "PF4")

ht <- sc_heatmap(
  pbmc,
  features = markers,
  label.features = label_genes,
  split.by = "cell_type",
  feature.groups = marker_groups,
  downsample = 30,
  seed = 2026,
  annotations = "cell_type",
  colors = sc_heatmap_palette("rdbu"),
  show.group.names = FALSE
)

ComplexHeatmap::draw(ht, merge_legend = TRUE)
save_sc_heatmap(ht, "pbmc3k_heatmap", width = 7, height = 5,
                merge_legend = TRUE)
```

如果需要更紧凑的组间比较，可以按一个或多个 metadata 字段计算平均表达：

```r
average_ht <- sc_heatmap(
  pbmc,
  features = markers,
  mode = "average",
  split.by = "cell_type",
  aggregate.by = "cell_type",
  feature.groups = marker_groups,
  annotations = "cell_type",
  show.column.names = TRUE
)
```

`order.by` 支持 pseudotime 等 metadata、marker 基因集合、自定义连续分数以及
`"cluster"`；原始文献中每个分组采用不同 marker 集合排序的方式继续通过
`sort.by` 支持。

`label.features` 会保留全部基因表达行，并以连线标注选定基因在热图中的位置。当
`annotations` 包含 `cell_type` 时，默认使用原 `divergentcolor` 配色；可通过
`annotation.colors` 显式覆盖。

连线基因标签默认使用斜体。设置 `label.fontface = "plain"` 可改为正体；也可用
`label.gp = grid::gpar(fontface = "plain", col = "red")` 完整控制样式。

不同显示层可以独立关闭：`show.feature.names` 控制基因标签，
`show.column.names` 控制细胞或聚合列名，`show.group.names` 控制分组切片标题，
`show.group.annotation` 控制顶部 metadata 注释。
每个热图分块默认带黑色边框。使用 `border.color = FALSE` 可关闭边框，也可以
传入其他颜色进行覆盖。

默认从当前 assay 的 `data` layer 读取表达值，逐基因进行 z-score，并截断到
`[-2, 2]`。在 `mode = "average"` 时，函数先在每个 `aggregate.by` 分组内，
对每个基因跨细胞计算算术平均值，再对聚合后的矩阵执行相同的标准化和截断。
设置 `scale = "none", clip = NULL` 可以查看未经标准化的分组均值。

在[使用案例集](inst/examples/scHeatmap-gallery.md)中可以查看五种完整使用方式。

所有参数请查看 `?sc_heatmap` 和 `?save_sc_heatmap`。

## 文献图复现

README 顶部的图片展示了本包最初针对的复杂场景：基因按生物学程序分组，同时
每个处理组使用不同的 marker 程序对组内细胞排序。复现数据来自 Lau 等发表于
Cell Reports 的 2020 年研究
[《IL-33-PU.1 Transcriptome Reprogramming Drives Functional State Transition
and Clearance Activity of Microglia in Alzheimer’s Disease》](https://www.cell.com/cell-reports/fulltext/S2211-1247(20)30430-7)
（Cell Reports 31，107530；DOI：10.1016/j.celrep.2020.107530）。通过 Git LFS
克隆完整仓库后，可以使用 `data-raw/` 中的 RDS 复现：

```r
cr_seu <- readRDS("data-raw/CR_2020_scRNAseq_seu_annoted.rds")
group_names <- c(WT = "WT", AD = "APP", IL33 = "APP+IL33")
cr_seu$treatment <- factor(
  unname(group_names[as.character(cr_seu$Group)]),
  levels = c("WT", "APP", "APP+IL33")
)
Seurat::Idents(cr_seu) <- "cell_type"
cr_subset <- subset(cr_seu, cell_type %in% c("Homeo", "DAM"))
SeuratObject::DefaultAssay(cr_subset) <- "RNA"

marker_groups <- list(
  `MHC-II` = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1"),
  Homeo = c("Cap1", "Fcrls", "Sft2d1", "Zfp69"),
  DAM = c("Nrp1", "Lpl", "Itgax", "Apoe", "Cst7")
)

cr_ht <- sc_heatmap(
  cr_subset,
  features = unlist(marker_groups, use.names = FALSE),
  split.by = "treatment",
  feature.groups = marker_groups,
  sort.by = list(
    WT = marker_groups[["MHC-II"]],
    APP = marker_groups$DAM,
    `APP+IL33` = marker_groups$Homeo
  ),
  decreasing = c(WT = FALSE, APP = FALSE, `APP+IL33` = TRUE),
  colors = sc_heatmap_palette("purple_black_yellow"),
  color.breaks = c(-1, 0, 1),
  clip = c(-1, 1),
  row_title_rot = 0,
  border = FALSE,
  row_names_gp = grid::gpar(fontface = "italic")
)

ComplexHeatmap::draw(cr_ht)
```

对应的图片导出代码以及 PBMC3K marker 自动筛选流程位于
[`development/demo.R`](development/demo.R)。

对于设置了固定热图主体尺寸的对象，可以计算完整图形所需的设备尺寸：

```r
size <- calc_ht_size(ht, unit = "inch")
save_sc_heatmap(ht, "my_heatmap", width = size["width"], height = size["height"])
```

## 内置配色

原项目 `utils.R` 中适合热图使用的配色已经整理到统一、带文档的接口中：

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

完整分析对象通过 Git LFS 保存在 `data-raw/`，不会被打包进安装后的 R 包。
原始复现分析和旧版工具函数保留在 `development/`，用于追溯和后续筛选。

## 开发与检查

以下命令应在仓库根目录运行。`load_all()` 会直接加载当前源码而不执行安装，适合
反复开发和调试：

```r
devtools::load_all(".")
devtools::document()
devtools::test()
devtools::check(document = FALSE, manual = FALSE, cran = FALSE)
```

在仓库根目录运行完整开发测试脚本：

```bash
conda run -n scptools Rscript development/demo.R
```

该脚本包含合成数据、PBMC3K marker 筛选以及 CR2020 文献图复现。完整运行需要
Seurat、dplyr、plyr、`pbmc3k.SeuratData`，以及通过 Git LFS 获取的
`data-raw/` RDS。如果只具备部分数据，也可以在交互式 R 会话中分段执行。

如果需要把当前源码正式安装到 `scptools` conda 环境的 R library，而不是仅使用
`load_all()` 临时加载，请运行：

```bash
conda run -n scptools Rscript -e 'devtools::install(".", upgrade = FALSE)'
conda run -n scptools Rscript -e 'library(scHeatmap); packageVersion("scHeatmap")'
```

这里的 `upgrade = FALSE` 只是不让 `devtools` 自动升级已经安装的 Seurat、
ComplexHeatmap 等依赖，不会阻止重新安装当前本地版本的 `scHeatmap`。由于 conda
环境中已经准备好了依赖，建议保留该设置，以免 R 在安装过程中意外升级依赖并
引入版本冲突。只有在明确希望安装时一并更新依赖的情况下，才使用
`upgrade = TRUE`。

安装后，只要激活 `scptools` 环境，就可以在任意工作目录使用
`library(scHeatmap)`。如需确认实际安装路径，可以查看该环境的 `.libPaths()`：

```bash
conda run -n scptools Rscript -e '.libPaths()'
```

## 许可证

MIT License。
