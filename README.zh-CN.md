# scHeatmap

[![English](https://img.shields.io/badge/lang-English-blue.svg)](README.md)
[![简体中文](https://img.shields.io/badge/语言-简体中文-red.svg)](README.zh-CN.md)

`scHeatmap` 是一个从 Seurat 对象直接创建单细胞热图的 R 包。它在
`ComplexHeatmap` 基础上提供简洁接口，同时明确保留基因顺序、marker 分组、
元数据分组以及组内细胞排序，使绘图过程更容易复现和调整。

![单细胞热图示例](man/figures/CR_2020_GSEreplicate_figS2E_20260721.jpg)

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

如果需要更紧凑的组间比较，可以按一个或多个 metadata 字段计算平均表达：

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
