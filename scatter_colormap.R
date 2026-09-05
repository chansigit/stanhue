# ══════════════════════════════════════════════════════════════════════════════
# scatter_colormap.R — 通用单细胞 UMAP 层级化自动配色工具 (R 版)
# ══════════════════════════════════════════════════════════════════════════════
#
# 根据任意 2D 降维坐标（UMAP / tSNE / PCA）和 cell type 标签，
# 自动推断 lineage 层级关系，生成"组间色相区分、组内明度渐变"的配色方案。
#
# 适用于任意物种、任意组织的单细胞数据。与 Python 版算法完全一致。
#
# 算法:
#   1. 计算每个 cell type 的 2D 质心
#   2. Ward 层次聚类 → 自动切分大类 (lineage)
#   3. 组内最大亚群排首位 (dominant)，其余按 dendrogram 叶序
#   4. 各大类以步长 2 在调色板上错开起始位 → dominant color 互不相同
#   5. 组内从偏移位顺延取色
#   6. 重叠感知 (overlap_aware=TRUE): 在 2D 空间中混在一起的类别，在各自
#      可选的槽位里挑 CIELAB ΔE 最远的颜色；类别数 ≤ 调色板长度时，
#      按最远点采样顺序发色（2 类不再是浅蓝 + 深蓝）
#
# Quick Start:
#   source("scatter_colormap.R")
#   colors <- assign_celltype_colors(umap_coords, cell_type_labels)
#   plot_umap(umap_coords, cell_type_labels, colors)
#
# 依赖: stats (base R), ggplot2, ggrepel (仅画图时需要)
# ══════════════════════════════════════════════════════════════════════════════


# ── 调色板 ──────────────────────────────────────────────────────────────────

#' ColorBrewer "Paired" 12 色调色板 — 6 对浅/深配色
#' 默认调色板，可通过 palette 参数替换。
PAIRED_PALETTE <- c(
  "#a6cee3",  #  1  浅蓝
  "#1f78b4",  #  2  深蓝
  "#b2df8a",  #  3  浅绿
  "#33a02c",  #  4  深绿
  "#fb9a99",  #  5  浅红
  "#e31a1c",  #  6  深红
  "#fdbf6f",  #  7  浅橙
  "#ff7f00",  #  8  深橙
  "#cab2d6",  #  9  浅紫
  "#6a3d9a",  # 10  深紫
  "#ffff99",  # 11  黄
  "#b15928"   # 12  棕
)


# ══════════════════════════════════════════════════════════════════════════════
#  输入验证
# ══════════════════════════════════════════════════════════════════════════════

#' 验证 coords 和 labels 的合法性
#' @param coords matrix, n x 2
#' @param labels character vector, length n
#' @return list(coords, labels) 规范化后的输入
.validate_inputs <- function(coords, labels) {
  coords <- as.matrix(coords)
  labels <- as.character(labels)
  errors <- character(0)

  if (ncol(coords) != 2)
    errors <- c(errors, sprintf("coords must have 2 columns, got %d", ncol(coords)))
  if (any(!is.finite(coords)))
    errors <- c(errors, "coords contains NaN, NA, or Inf values")
  if (length(labels) != nrow(coords))
    errors <- c(errors, sprintf("labels length (%d) != coords rows (%d)",
                                length(labels), nrow(coords)))
  if (length(errors) > 0)
    stop("Input validation failed:\n",
         paste("  -", errors, collapse = "\n"), call. = FALSE)

  list(coords = coords, labels = labels)
}


# ══════════════════════════════════════════════════════════════════════════════
#  核心 API
# ══════════════════════════════════════════════════════════════════════════════

#' 共享聚类流水线（内部函数）
#' @return list(groups, unique_types, type_counts, coords, codes)
.cluster_pipeline <- function(coords, labels, n_major_groups = NULL) {
  validated <- .validate_inputs(coords, labels)
  coords <- validated$coords
  labels <- validated$labels

  unique_types <- sort(unique(labels))
  n_types <- length(unique_types)
  codes <- as.integer(factor(labels, levels = unique_types))  # 1-indexed

  if (n_types < 2) {
    groups <- list()
    type_counts <- integer(0)
    if (n_types == 1) {
      groups[["1"]] <- unique_types
      type_counts <- setNames(sum(labels == unique_types), unique_types)
    }
    return(list(groups = groups, unique_types = unique_types,
                type_counts = type_counts, coords = coords, codes = codes))
  }

  # 向量化质心计算
  centroids <- .compute_centroids_fast(coords, labels, unique_types)
  hc <- hclust(dist(centroids), method = "ward.D2")

  if (!is.null(n_major_groups) && n_major_groups < 1)
    stop("n_major_groups must be >= 1, got ", n_major_groups, call. = FALSE)
  if (is.null(n_major_groups)) {
    n_major_groups <- .auto_determine_k(hc, n_types)
  }
  n_major_groups <- min(n_major_groups, n_types)

  cluster_ids <- cutree(hc, k = n_major_groups)
  type_counts <- table(labels)[unique_types]
  groups <- .build_ordered_groups(unique_types, cluster_ids, type_counts, hc)
  list(groups = groups, unique_types = unique_types,
       type_counts = type_counts, coords = coords, codes = codes)
}


#' 给定 2D 坐标和细胞标签，返回 named vector: cell_type → hex_color
#'
#' @param coords     matrix/data.frame, n_cells x 2, 任意 2D 降维坐标
#' @param labels     character vector, 长度 n_cells, 每个细胞的类型标签
#' @param n_major_groups  integer|NULL, 手动指定大类数量。NULL 自动确定
#' @param palette    character vector, 有序调色板。NULL 使用 PAIRED_PALETTE
#' @param return_groups  logical, 若 TRUE 返回 list(colors, groups)
#' @param overlap_aware  logical, 默认 TRUE。在 2D 空间中相互重叠的类别会被
#'                   分配感知距离 (CIELAB ΔE) 尽量远的颜色；不重叠的类别
#'                   保持原有层级配色。FALSE 完全回退到旧行为
#' @param overlap_threshold  numeric [0,1], 默认 0.1。两类别的空间混合度
#'                   (见 .compute_overlap) 低于该值时视为"不重叠"
#' @return named character vector: cell_type → hex color;
#'         若 return_groups=TRUE, 返回 list(colors=..., groups=...)
assign_celltype_colors <- function(coords, labels,
                                   n_major_groups = NULL,
                                   palette = NULL,
                                   return_groups = FALSE,
                                   overlap_aware = TRUE,
                                   overlap_threshold = 0.1) {
  if (is.null(palette)) palette <- PAIRED_PALETTE
  n_pal <- length(palette)

  res <- .cluster_pipeline(coords, labels, n_major_groups)
  groups <- res$groups
  unique_types <- res$unique_types
  type_counts <- res$type_counts
  n_types <- length(unique_types)

  # 边界情况
  .ret <- function(cm) {
    if (return_groups) list(colors = cm, groups = groups) else cm
  }
  if (n_types == 0) return(.ret(character(0)))
  if (n_types == 1) return(.ret(setNames(palette[1], unique_types)))

  # 重叠矩阵 + 调色板内两两感知色差
  if (overlap_aware) {
    overlap <- .compute_overlap(res$coords, res$codes, n_types)
    overlap[overlap < overlap_threshold] <- 0
    delta_e <- .palette_delta_e(palette)
  } else {
    overlap <- matrix(0, n_types, n_types)
    delta_e <- matrix(0, n_pal, n_pal)
  }
  type_to_idx <- setNames(seq_len(n_types), unique_types)
  assigned <- rep(NA_integer_, n_types)  # type_idx → palette_idx (1-indexed)

  if (n_types <= n_pal && is.null(n_major_groups)) {
    # 类型数 ≤ 调色板长度：每个类型独占一色（保证零重复）。
    # overlap_aware 时按 CIELAB 最远点采样顺序发色，先发差异最大的颜色
    slots <- if (overlap_aware) .spread_order(palette, delta_e) else seq_len(n_types)
    assigned <- .greedy_assign(seq_len(n_types), slots, overlap, delta_e, assigned)
    return(.ret(setNames(palette[assigned], unique_types)))
  }

  # Step 4: 按总细胞数排序，分配偏移
  group_totals <- vapply(groups, function(members) {
    sum(type_counts[members])
  }, numeric(1))
  sorted_gids <- names(sort(group_totals, decreasing = TRUE))
  n_groups <- length(sorted_gids)

  if (n_groups <= n_pal %/% 2) {
    offsets <- seq(0, by = 2, length.out = n_groups)   # 步长 2
  } else {
    base <- c(seq(0, n_pal - 1, by = 2), seq(1, n_pal - 1, by = 2))
    # cycle when n_groups > n_pal
    offsets <- base[((seq_len(n_groups) - 1) %% length(base)) + 1]
  }

  # Step 5: 取色，两级贪心（无重叠时与旧算法完全一致）：
  #   组级：每组 dominant 从剩余 offset 池中挑与已分配重叠类别色差最大的位置
  #        （无重叠则保留自己原本的 offset），全组占据从该 offset 起的连续一段 → 色系结构不变
  #   组内：其余成员若与已分配类别重叠，则在本组槽位中挑感知距离最远的，
  #        否则按叶序顺延
  # offset 池用 1-indexed palette 位置表示
  dominants <- vapply(sorted_gids, function(gid) unname(type_to_idx[groups[[gid]][1]]), integer(1))
  assigned <- .greedy_assign(unname(dominants), offsets + 1, overlap, delta_e, assigned)
  for (gid in sorted_gids) {
    members <- unname(type_to_idx[groups[[gid]]])
    start <- assigned[members[1]]
    slots <- ((start - 1 + seq_along(members) - 1) %% n_pal) + 1
    assigned <- .greedy_assign(members[-1], slots[-1], overlap, delta_e, assigned)
  }

  .ret(setNames(palette[assigned], unique_types))
}


#' 只返回分组结果（不分配颜色），用于检查 / 调试
#'
#' @inheritParams assign_celltype_colors
#' @return named list: group_id → character vector of cell types (首位是 dominant)
get_groups <- function(coords, labels, n_major_groups = NULL) {
  .cluster_pipeline(coords, labels, n_major_groups)$groups
}


# ══════════════════════════════════════════════════════════════════════════════
#  可视化
# ══════════════════════════════════════════════════════════════════════════════

#' 绘制出版级 UMAP scatter plot (ggplot2)
#'
#' @param coords     matrix, n_cells x 2
#' @param labels     character vector, n_cells
#' @param color_map  named character vector (assign_celltype_colors 的输出)
#' @param groups     list|NULL, 分组 (get_groups 的输出)。NULL 自动推断
#' @param point_size numeric, 点大小
#' @param show_group_labels logical, 是否标注大群名
#' @param label_size numeric, 标注字号
#' @param xlabel,ylabel character, 坐标轴标签
#' @param title      character|NULL, 标题
#' @param legend_ncol integer|NULL, 图例列数。NULL 自动
#' @param save_path  character|NULL, 保存路径
#' @param width,height,dpi  保存参数
#' @return ggplot 对象 (invisible)
plot_umap <- function(coords, labels, color_map,
                      groups = NULL,
                      point_size = 0.3,
                      show_group_labels = TRUE,
                      label_size = 4,
                      xlabel = "UMAP 1",
                      ylabel = "UMAP 2",
                      title = NULL,
                      legend_ncol = NULL,
                      save_path = NULL,
                      width = 14, height = 9, dpi = 300) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot_umap(). Install with: install.packages('ggplot2')")

  coords <- as.matrix(coords)
  labels <- as.character(labels)
  cell_counts <- table(labels)

  if (is.null(groups)) groups <- get_groups(coords, labels)

  # 按总细胞数排序大群
  group_totals <- vapply(groups, function(m) sum(cell_counts[m]), numeric(1))
  sorted_gids <- names(sort(group_totals, decreasing = TRUE))
  group_names <- vapply(groups[sorted_gids], .group_display_name, character(1))
  names(group_names) <- sorted_gids

  # 构建 data.frame，按细胞数排序（稀有群画在上面）
  df <- data.frame(
    x     = coords[, 1],
    y     = coords[, 2],
    ct    = labels,
    stringsAsFactors = FALSE
  )

  # 建立 factor level 顺序：大群先画 → 稀有群后画（ggplot 后画的在上层）
  ct_order <- unlist(lapply(rev(sorted_gids), function(gid) {
    members <- groups[[gid]]
    members[order(cell_counts[members], decreasing = TRUE)]
  }))
  # 反转：让 count 小的排后面（后画）
  ct_order_plot <- rev(ct_order)
  df$ct <- factor(df$ct, levels = ct_order_plot)

  # 配色 vector (与 factor levels 对齐)
  pal <- color_map[levels(df$ct)]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = ct)) +
    ggplot2::geom_point(size = point_size, stroke = 0, shape = 16) +
    ggplot2::scale_color_manual(values = pal, name = NULL) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = xlabel, y = ylabel, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      legend.key.size  = ggplot2::unit(3, "mm"),
      legend.text      = ggplot2::element_text(size = 7),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3),
        ncol = if (!is.null(legend_ncol)) legend_ncol
               else max(1, ceiling(length(color_map) / 20))
      )
    )

  # 质心标注
  if (show_group_labels) {
    label_df <- do.call(rbind, lapply(sorted_gids, function(gid) {
      mask <- labels %in% groups[[gid]]
      data.frame(
        x    = mean(coords[mask, 1]),
        y    = mean(coords[mask, 2]),
        name = group_names[gid],
        stringsAsFactors = FALSE
      )
    }))

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = name),
        inherit.aes = FALSE,
        size = label_size, fontface = "bold", color = "#222222",
        bg.color = "white", bg.r = 0.15,
        point.size = NA, box.padding = 0.3,
        max.overlaps = Inf, seed = 42
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = name),
        inherit.aes = FALSE,
        size = label_size, fontface = "bold", color = "#222222"
      )
    }
  }

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height, dpi = dpi,
                    bg = "white")
    message("Saved -> ", save_path)
  }

  print(p)
  invisible(p)
}


#' 打印分组配色摘要到 console
#'
#' @param color_map  named character vector
#' @param groups     list (get_groups 的输出)
#' @param labels     character vector (原始标签，用于统计细胞数)
print_color_summary <- function(color_map, groups, labels) {
  labels <- as.character(labels)
  cell_counts <- table(labels)

  group_totals <- vapply(groups, function(m) sum(cell_counts[m]), numeric(1))
  sorted_gids <- names(sort(group_totals, decreasing = TRUE))

  n_cells <- length(labels)
  n_types <- length(color_map)
  n_groups <- length(groups)
  cat(sprintf("\n  %s cells | %d cell types | %d groups\n",
              format(n_cells, big.mark = ","), n_types, n_groups))
  cat(strrep("=", 60), "\n")

  for (gid in sorted_gids) {
    members <- groups[[gid]]
    total <- sum(cell_counts[members])
    gname <- .group_display_name(members)
    cat(sprintf("\n  -- %s (%d subtypes, %s cells) --\n",
                gname, length(members), format(total, big.mark = ",")))
    for (i in seq_along(members)) {
      ct <- members[i]
      star <- if (i == 1) "\u2605" else " "
      cat(sprintf("   %s %-28s  %s  (%s)\n",
                  star, ct, color_map[ct],
                  format(cell_counts[ct], big.mark = ",")))
    }
  }
  cat("\n", strrep("=", 60), "\n")
  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
#  Seurat / SingleCellExperiment 便捷接口
# ══════════════════════════════════════════════════════════════════════════════

#' 从 Seurat 对象自动提取坐标和标签并配色
#'
#' @param seurat_obj Seurat 对象
#' @param reduction  character, 降维名称 (默认 "umap")
#' @param group_by   character, metadata 列名 (默认 "cell_type")
#' @param ...        传给 assign_celltype_colors 的其他参数
#' @return named character vector: cell_type → hex color
color_seurat <- function(seurat_obj, reduction = "umap",
                         group_by = "cell_type", ...) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Seurat is required. Install with: install.packages('Seurat')")

  coords <- Seurat::Embeddings(seurat_obj, reduction = reduction)[, 1:2]
  labels <- seurat_obj[[group_by, drop = TRUE]]
  assign_celltype_colors(coords, labels, ...)
}


#' 从 SingleCellExperiment 对象自动提取坐标和标签并配色
#'
#' @param sce        SingleCellExperiment 对象
#' @param dimred     character, 降维名称 (默认 "UMAP")
#' @param col_name   character, colData 列名 (默认 "cell_type")
#' @param ...        传给 assign_celltype_colors 的其他参数
#' @return named character vector: cell_type → hex color
color_sce <- function(sce, dimred = "UMAP", col_name = "cell_type", ...) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("SingleCellExperiment is required.")

  coords <- SingleCellExperiment::reducedDim(sce, dimred)[, 1:2]
  labels <- SummarizedExperiment::colData(sce)[[col_name]]
  assign_celltype_colors(coords, labels, ...)
}


# ══════════════════════════════════════════════════════════════════════════════
#  内部工具函数
# ══════════════════════════════════════════════════════════════════════════════

.compute_centroids_fast <- function(coords, labels, unique_types) {
  # 向量化：用 factor + rowsum 代替逐类 Python 循环
  f <- factor(labels, levels = unique_types)
  codes <- as.integer(f)
  sums <- rowsum(coords, codes, reorder = FALSE)
  counts <- tabulate(codes, nbins = length(unique_types))
  counts[counts == 0] <- 1  # 避免除零
  centroids <- sums / counts
  rownames(centroids) <- unique_types
  centroids
}


.compute_overlap <- function(coords, codes, n_types, n_bins = NULL) {
  # 类别两两之间的 2D 空间重叠度矩阵，取值 [0, 1]，对角线为 0。
  # 把坐标切成 n_bins x n_bins 网格。对类别 a 的每个细胞，看它所在格子里
  # 属于类别 b 的细胞占比，再对 a 的所有细胞取平均 —— 网格近似的
  # "近邻标签混合度" mix[a, b]。取 max(mix[a, b], mix[b, a]) 对称化，
  # 这样"小群嵌在大群里"也能被识别。
  # 网格粗细随细胞数自适应：平均每格约 20 个细胞，范围 [16, 128]。
  n_cells <- nrow(coords)
  if (is.null(n_bins))
    n_bins <- as.integer(min(128, max(16, round(sqrt(n_cells / 20)))))

  mins <- apply(coords, 2, min)
  spans <- apply(coords, 2, max) - mins
  spans[spans == 0] <- 1
  ix <- pmin(floor((coords[, 1] - mins[1]) / spans[1] * n_bins), n_bins - 1)
  iy <- pmin(floor((coords[, 2] - mins[2]) / spans[2] * n_bins), n_bins - 1)
  bin_id <- ix * n_bins + iy + 1  # 1-indexed
  n_total_bins <- n_bins * n_bins

  # H[t, bin] = 类别 t 在该格子的细胞数
  H <- matrix(
    tabulate((codes - 1) * n_total_bins + bin_id, nbins = n_types * n_total_bins),
    nrow = n_types, byrow = TRUE
  )
  share <- sweep(H, 2, pmax(colSums(H), 1), "/")  # 每格内各类别占比
  p <- H / pmax(rowSums(H), 1)                    # 各类别的空间分布
  mix <- p %*% t(share)                           # mix[a, b]
  overlap <- pmax(mix, t(mix))
  diag(overlap) <- 0
  overlap
}


.hex_to_lab <- function(hex_colors) {
  # hex (sRGB, D65) → CIELAB，n x 3 矩阵。与 Python 版公式一致。
  rgb <- t(grDevices::col2rgb(hex_colors)) / 255
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  m <- matrix(c(0.4124564, 0.3575761, 0.1804375,
                0.2126729, 0.7151522, 0.0721750,
                0.0193339, 0.1191920, 0.9503041), nrow = 3, byrow = TRUE)
  xyz <- sweep(lin %*% t(m), 2, c(0.95047, 1.0, 1.08883), "/")
  d <- 6 / 29
  f <- ifelse(xyz > d^3, xyz^(1 / 3), xyz / (3 * d^2) + 4 / 29)
  cbind(L = 116 * f[, 2] - 16,
        a = 500 * (f[, 1] - f[, 2]),
        b = 200 * (f[, 2] - f[, 3]))
}


.palette_delta_e <- function(palette) {
  # 调色板内两两 CIE76 ΔE 矩阵
  as.matrix(dist(.hex_to_lab(palette)))
}


.spread_order <- function(palette, delta_e) {
  # 调色板的"最远点采样"顺序：以白色 (L=100, a=b=0) 为第 0 个锚点，
  # 第一个选中的是与白色对比最强的深饱和色，之后每次选与已选颜色最小 ΔE 最大的。
  n <- nrow(delta_e)
  lab <- .hex_to_lab(palette)
  min_d <- sqrt(rowSums(sweep(lab, 2, c(100, 0, 0))^2))
  start <- which.max(min_d)
  ord <- start
  chosen <- logical(n)
  chosen[start] <- TRUE
  min_d <- pmin(min_d, delta_e[start, ])
  for (k in seq_len(n - 1)) {
    min_d[chosen] <- -Inf
    nxt <- which.max(min_d)
    ord <- c(ord, nxt)
    chosen[nxt] <- TRUE
    min_d <- pmin(min_d, delta_e[nxt, ])
  }
  unname(ord)
}


.greedy_assign <- function(type_ids, slots, overlap, delta_e, assigned) {
  # 重叠感知的贪心取色。assigned: 长度 n_types 的整数向量 (NA = 未分配)，
  # 返回更新后的向量。
  # - slots 是候选调色板位置（按默认优先顺序），每个位置只用一次
  # - 处理顺序：与其他类别总重叠量大的先选
  # - 类别 a 在位置 s 的得分 = sum_b overlap[a, b] * ΔE(s, color_b)，
  #   b 遍历已分配类别。得分全 0（无重叠）→ 取自己的默认位置 slots[i]
  #   （若已被占则取第一个未用位置）。
  if (length(type_ids) == 0) return(assigned)
  if (length(slots) < length(type_ids))
    stop(sprintf("need %d slots but only %d given", length(type_ids), length(slots)))
  remaining <- slots
  # 默认位置：type_ids[i] ↔ slots[i]（旧算法的结果）。无重叠的类别尽量保留
  # 默认位置，避免一个类别换位后其余类别全部顺延（级联）。
  preferred <- setNames(slots[seq_along(type_ids)], type_ids)
  total_overlap <- rowSums(overlap)
  ord <- type_ids[order(-total_overlap[type_ids], seq_along(type_ids))]

  for (a in ord) {
    pref <- preferred[[as.character(a)]]
    best <- if (pref %in% remaining) pref else remaining[1]
    done <- which(!is.na(assigned))
    if (length(done) > 0) {
      w <- overlap[a, done]
      if (any(w > 0)) {
        scores <- delta_e[remaining, assigned[done], drop = FALSE] %*% w
        best <- remaining[which.max(scores)]
      }
    }
    assigned[a] <- best
    remaining <- remaining[-match(best, remaining)]
  }
  assigned
}


.auto_determine_k <- function(hc, n_types) {
  # 在 dendrogram 合并距离中寻找自然切分点
  # 策略：找第一个超过中位数 2 倍的显著跳跃，避免被末尾巨大 gap 吸引
  if (n_types <= 3) return(n_types)

  heights <- hc$height
  gaps <- diff(heights)
  rel_gaps <- gaps / (heights[-length(heights)] + 1e-10)

  k_max <- max(15, n_types %/% 4)
  k_min <- 3
  n_gaps <- length(rel_gaps)

  # gap index i (1-based) 对应 k = n_types - i
  idx_hi <- n_gaps - k_min + 1       # k=k_min
  idx_lo <- max(1, n_gaps - k_max + 1) # k=k_max
  if (idx_lo > idx_hi) { tmp <- idx_lo; idx_lo <- idx_hi; idx_hi <- tmp }

  window_vals <- rel_gaps[idx_lo:idx_hi]
  if (length(window_vals) == 0) return(k_min)

  median_gap <- median(window_vals)
  threshold <- median_gap * 2.0

  # 回退值：窗口内全局最大
  best_idx <- idx_lo - 1 + which.max(window_vals)

  # 从低 index (大 k) 往高 index (小 k) 扫描
  for (i in idx_lo:idx_hi) {
    if (rel_gaps[i] >= threshold) {
      best_idx <- i
      break
    }
  }

  k <- n_types - best_idx
  max(k_min, min(k, k_max))
}


.build_ordered_groups <- function(unique_types, cluster_ids, type_counts, hc) {
  # cluster_ids: named vector; type_counts: named vector aligned with unique_types
  groups <- split(unique_types, cluster_ids)

  # dendrogram 叶序
  leaf_order <- hc$order  # 1-indexed positions
  type_order <- setNames(seq_along(leaf_order), unique_types[leaf_order])

  groups <- lapply(groups, function(members) {
    # dominant = 细胞数最多的
    counts_m <- type_counts[members]
    dominant <- names(which.max(counts_m))
    rest <- setdiff(members, dominant)
    rest <- rest[order(type_order[rest])]
    c(dominant, rest)
  })

  groups
}


.common_prefix <- function(strings) {
  if (length(strings) == 0) return("")
  prefix <- strings[1]
  for (s in strings[-1]) {
    while (!startsWith(s, prefix)) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
      if (nchar(prefix) == 0) return("")
    }
  }
  # 去掉尾部空格/下划线
  clean <- sub("[_ ]+$", "", prefix)
  if (nchar(clean) == 0) return("")

  # 如果末尾是数字，检查是否截在数字中间
  last_char <- substr(clean, nchar(clean), nchar(clean))
  if (grepl("[0-9]", last_char)) {
    for (s in strings) {
      if (nchar(s) > nchar(clean)) {
        next_char <- substr(s, nchar(clean) + 1, nchar(clean) + 1)
        if (grepl("[0-9]", next_char)) {
          # 截在数字中间，回退到最后一个非数字字符
          m <- regexpr("[^0-9][0-9]+$", clean)
          if (m > 0) {
            candidate <- sub("[_ ]+$", "", substr(clean, 1, m))
            return(if (nchar(candidate) >= 2) candidate else "")
          }
          return("")
        }
      }
    }
  }
  clean
}


.group_display_name <- function(members) {
  # 公共前缀 → 多数票首词 → dominant 名
  prefix <- .common_prefix(members)
  if (nchar(prefix) >= 3) return(prefix)

  # 多数票首词
  first_tokens <- sub("[_ ].*$", "", members)
  if (length(first_tokens) > 0) {
    tt <- sort(table(first_tokens), decreasing = TRUE)
    if (tt[1] > length(members) * 0.4) {
      return(names(tt)[1])
    }
  }
  members[1]
}
