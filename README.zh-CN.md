# scHeatmap

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![简体中文](https://img.shields.io/badge/语言-简体中文-red.svg)](README.zh-CN.md)

`scHeatmap` 是一个从 Seurat 对象直接创建单细胞热图的 R 包。它在
`ComplexHeatmap` 基础上提供简洁接口，同时明确保留基因顺序、marker 分组、
元数据分组以及组内细胞排序，使绘图过程更容易复现和调整。

![单细胞热图示例](man/figures/CR_2020_GSEreplicate_figS2E_20260721.jpg)

## 主要功能

- 支持读取 Seurat v5 assay 和 layer。
- 使用任意 Seurat metadata 字段对细胞分组。
- 按指定顺序组织基因和 marker 集合。
- 根据 marker 集合得分对每组细胞分别排序。
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

默认情况下，函数从当前 assay 的 `data` layer 读取表达值，对每个基因在所选
细胞中进行 z-score 标准化，并将显示范围截断至 `[-2, 2]`。如果希望直接展示
指定 layer 的数值，可以设置 `scale = FALSE, clip = NULL`。

所有参数请查看 `?sc_heatmap` 和 `?save_sc_heatmap`。

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

```r
devtools::document()
devtools::test()
devtools::check()
```

## 许可证

MIT License。
