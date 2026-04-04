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

#' 给定 2D 坐标和细胞标签，返回 named vector: cell_type → hex_color
#'
#' @param coords     matrix/data.frame, n_cells x 2, 任意 2D 降维坐标
#' @param labels     character vector, 长度 n_cells, 每个细胞的类型标签
#' @param n_major_groups  integer|NULL, 手动指定大类数量。NULL 自动确定
#' @param palette    character vector, 有序调色板。NULL 使用 PAIRED_PALETTE
#' @return named character vector: cell_type → hex color
assign_celltype_colors <- function(coords, labels,
                                   n_major_groups = NULL,
                                   palette = NULL) {
  validated <- .validate_inputs(coords, labels)
  coords <- validated$coords
  labels <- validated$labels

  if (is.null(palette)) palette <- PAIRED_PALETTE
  n_pal <- length(palette)

  unique_types <- sort(unique(labels))
  n_types <- length(unique_types)

  # 边界情况
  if (n_types == 0) return(character(0))
  if (n_types == 1) return(setNames(palette[1], unique_types))
  if (n_types <= n_pal && is.null(n_major_groups)) {
    return(setNames(palette[seq_len(n_types)], unique_types))
  }

  # Step 1-3: 聚类 + 分组 + 组内排序
  centroids <- .compute_centroids(coords, labels, unique_types)
  hc <- hclust(dist(centroids), method = "ward.D2")

  if (!is.null(n_major_groups) && n_major_groups < 1)
    stop("n_major_groups must be >= 1, got ", n_major_groups, call. = FALSE)
  if (is.null(n_major_groups)) {
    n_major_groups <- .auto_determine_k(hc, n_types)
  }
  n_major_groups <- min(n_major_groups, n_types)

  cluster_ids <- cutree(hc, k = n_major_groups)  # named by row
  groups <- .build_ordered_groups(unique_types, cluster_ids, labels, hc)

  # Step 4: 按总细胞数排序，分配偏移
  cell_counts <- table(labels)
  group_totals <- vapply(groups, function(members) {
    sum(cell_counts[members])
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

  # Step 5: 组内顺延取色 (0-indexed offset, 转 1-indexed palette)
  color_map <- character(0)
  for (idx in seq_along(sorted_gids)) {
    gid <- sorted_gids[idx]
    offset <- offsets[idx]
    members <- groups[[gid]]
    for (i in seq_along(members)) {
      pal_idx <- ((offset + i - 1) %% n_pal) + 1  # R 的 1-indexed
      color_map[members[i]] <- palette[pal_idx]
    }
  }

  return(color_map)
}


#' 只返回分组结果（不分配颜色），用于检查 / 调试
#'
#' @inheritParams assign_celltype_colors
#' @return named list: group_id → character vector of cell types (首位是 dominant)
get_groups <- function(coords, labels, n_major_groups = NULL) {
  validated <- .validate_inputs(coords, labels)
  coords <- validated$coords
  labels <- validated$labels
  unique_types <- sort(unique(labels))
  n_types <- length(unique_types)

  centroids <- .compute_centroids(coords, labels, unique_types)
  hc <- hclust(dist(centroids), method = "ward.D2")

  if (!is.null(n_major_groups) && n_major_groups < 1)
    stop("n_major_groups must be >= 1, got ", n_major_groups, call. = FALSE)
  if (is.null(n_major_groups)) {
    n_major_groups <- .auto_determine_k(hc, n_types)
  }
  n_major_groups <- min(n_major_groups, n_types)

  cluster_ids <- cutree(hc, k = n_major_groups)
  .build_ordered_groups(unique_types, cluster_ids, labels, hc)
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

.compute_centroids <- function(coords, labels, unique_types) {
  centroids <- matrix(0, nrow = length(unique_types), ncol = 2)
  rownames(centroids) <- unique_types
  for (i in seq_along(unique_types)) {
    mask <- labels == unique_types[i]
    centroids[i, ] <- colMeans(coords[mask, , drop = FALSE])
  }
  centroids
}


.auto_determine_k <- function(hc, n_types) {
  # 在 dendrogram 合并距离中找最大相对间距跳跃
  if (n_types <= 3) return(n_types)

  heights <- hc$height  # 已排序（升序）
  gaps <- diff(heights)
  rel_gaps <- gaps / (heights[-length(heights)] + 1e-10)

  n_gaps <- length(rel_gaps)
  search_start <- max(1, n_gaps - min(25, n_types - 1) + 1)
  best_idx <- search_start - 1 + which.max(rel_gaps[search_start:n_gaps])

  # best_idx 是在 heights 中的位置；对应切 k 个簇
  k <- n_types - best_idx
  max(3, min(k, 15))
}


.build_ordered_groups <- function(unique_types, cluster_ids, all_labels, hc) {
  # cluster_ids: named vector (names = unique_types 的行号对应名称)
  groups <- split(unique_types, cluster_ids)

  cell_counts <- table(all_labels)
  # dendrogram 叶序
  leaf_order <- hc$order  # 1-indexed positions
  type_order <- setNames(seq_along(leaf_order), unique_types[leaf_order])

  groups <- lapply(groups, function(members) {
    # dominant = 细胞数最多的
    counts_m <- cell_counts[members]
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
