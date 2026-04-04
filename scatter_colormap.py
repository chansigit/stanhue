"""
scatter_colormap.py — 通用单细胞 UMAP 层级化自动配色工具
=====================================================

根据任意 2D 降维坐标（UMAP / tSNE / PCA）和 cell type 标签，
自动推断 lineage 层级关系，生成"组间色相区分、组内明度渐变"的配色方案。

适用于任意物种、任意组织的单细胞数据。

算法:
  1. 计算每个 cell type 的 2D 质心
  2. Ward 层次聚类 → 自动切分大类 (lineage)
  3. 组内最大亚群排首位 (dominant)，其余按 dendrogram 叶序
  4. 各大类以步长 2 在调色板上错开起始位 → dominant color 互不相同
  5. 组内从偏移位顺延取色

Quick Start:
  >>> from scatter_colormap import assign_celltype_colors, plot_umap
  >>> colors = assign_celltype_colors(adata.obsm["X_umap"], adata.obs["cell_type"])
  >>> plot_umap(adata.obsm["X_umap"], adata.obs["cell_type"], colors, save_path="umap.png")

依赖: numpy, scipy, matplotlib (仅画图时需要)
"""

from __future__ import annotations

import re
import numpy as np
from collections import Counter
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list


# ══════════════════════════════════════════════════════════
#  调色板
# ══════════════════════════════════════════════════════════

#: ColorBrewer "Paired" 12 色调色板 — 6 对浅/深配色。
#: 默认调色板，可通过 ``palette`` 参数替换为任意列表。
PAIRED_PALETTE: list[str] = [
    "#a6cee3",  #  0  浅蓝
    "#1f78b4",  #  1  深蓝
    "#b2df8a",  #  2  浅绿
    "#33a02c",  #  3  深绿
    "#fb9a99",  #  4  浅红
    "#e31a1c",  #  5  深红
    "#fdbf6f",  #  6  浅橙
    "#ff7f00",  #  7  深橙
    "#cab2d6",  #  8  浅紫
    "#6a3d9a",  #  9  深紫
    "#ffff99",  # 10  黄
    "#b15928",  # 11  棕
]


# ══════════════════════════════════════════════════════════
#  输入验证
# ══════════════════════════════════════════════════════════

def _validate_inputs(
    coords: np.ndarray, labels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """验证输入合法性，返回规范化的 (coords, labels)。"""
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    errors: list[str] = []
    if coords.ndim != 2 or coords.shape[1] != 2:
        errors.append(f"coords must be shape (n, 2), got {coords.shape}")
    if coords.size > 0 and np.any(~np.isfinite(coords)):
        errors.append("coords contains NaN or Inf values")
    if coords.ndim == 2 and len(labels) != coords.shape[0]:
        errors.append(
            f"labels length ({len(labels)}) != coords rows ({coords.shape[0]})"
        )
    if errors:
        raise ValueError(
            "Input validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        )
    return coords, labels


# ══════════════════════════════════════════════════════════
#  核心 API
# ══════════════════════════════════════════════════════════

def _cluster_pipeline(
    coords: np.ndarray,
    labels: np.ndarray,
    n_major_groups: int | None = None,
) -> tuple[dict[int, list[str]], np.ndarray, np.ndarray, np.ndarray]:
    """共享聚类流水线，返回 (groups, unique_types, type_counts, Z)。"""
    coords, labels = _validate_inputs(coords, labels)

    # 向量化：字符串标签 → 整数编码，后续全部用整数操作
    unique_types, codes = np.unique(labels, return_inverse=True)
    n_types = len(unique_types)

    if n_types < 2:
        groups: dict[int, list[str]] = {}
        if n_types == 1:
            groups[1] = [unique_types[0]]
        type_counts = np.bincount(codes, minlength=n_types)
        return groups, unique_types, type_counts, np.empty((0, 4))

    # 向量化质心计算
    centroids = _compute_centroids_fast(coords, codes, n_types)
    Z = linkage(centroids, method="ward")

    if n_major_groups is not None and n_major_groups < 1:
        raise ValueError(f"n_major_groups must be >= 1, got {n_major_groups}")
    if n_major_groups is None:
        n_major_groups = _auto_determine_k(Z, n_types)
    n_major_groups = min(n_major_groups, n_types)

    major_labels = fcluster(Z, t=n_major_groups, criterion="maxclust")
    type_counts = np.bincount(codes, minlength=n_types)
    groups = _build_ordered_groups(unique_types, major_labels, type_counts, Z)
    return groups, unique_types, type_counts, Z


def assign_celltype_colors(
    coords: np.ndarray,
    labels: np.ndarray,
    n_major_groups: int | None = None,
    palette: list[str] | None = None,
    return_groups: bool = False,
) -> dict[str, str] | tuple[dict[str, str], dict[int, list[str]]]:
    """
    给定 2D 坐标和细胞标签，返回 ``{cell_type: hex_color}``。

    Parameters
    ----------
    coords : ndarray, shape (n_cells, 2)
        任意 2D 降维坐标（UMAP / tSNE / PCA 均可）。
    labels : array-like, shape (n_cells,)
        每个细胞的类型标签（字符串）。
    n_major_groups : int | None
        手动指定大类数量。``None`` 则由算法自动确定。
    palette : list[str] | None
        组内有序调色板。``None`` 使用 ``PAIRED_PALETTE``。
        可传入任意长度的 hex 颜色列表。
    return_groups : bool
        若为 True，同时返回分组信息，避免重复聚类。

    Returns
    -------
    dict[str, str]  或  (dict[str, str], dict[int, list[str]])
        cell_type → hex color 映射。若 return_groups=True，额外返回分组。
    """
    if palette is None:
        palette = PAIRED_PALETTE
    n_pal = len(palette)

    groups, unique_types, type_counts, Z = _cluster_pipeline(
        coords, labels, n_major_groups,
    )
    n_types = len(unique_types)

    # 边界情况
    if n_types == 0:
        return ({}, {}) if return_groups else {}
    if n_types == 1:
        cm = {unique_types[0]: palette[0]}
        return (cm, groups) if return_groups else cm
    if n_types <= n_pal and n_major_groups is None:
        # 类型数 ≤ 调色板长度：每个类型独占一色（保证零重复）。
        # groups 仍来自空间聚类，可用于图例分组 / 标签标注。
        cm = {ct: palette[i] for i, ct in enumerate(unique_types)}
        return (cm, groups) if return_groups else cm

    # Step 4: 按总细胞数排序，分配偏移
    sorted_gids = sorted(
        groups,
        key=lambda g: -sum(type_counts[np.searchsorted(unique_types, ct)]
                           for ct in groups[g]),
    )
    n_groups = len(sorted_gids)

    if n_groups <= n_pal // 2:
        offsets = [i * 2 for i in range(n_groups)]
    else:
        base = list(range(0, n_pal, 2)) + list(range(1, n_pal, 2))
        offsets = [base[i % len(base)] for i in range(n_groups)]

    # Step 5: 组内顺延取色
    color_map: dict[str, str] = {}
    for gid, offset in zip(sorted_gids, offsets):
        for i, ct in enumerate(groups[gid]):
            color_map[ct] = palette[(offset + i) % n_pal]

    return (color_map, groups) if return_groups else color_map


def get_groups(
    coords: np.ndarray,
    labels: np.ndarray,
    n_major_groups: int | None = None,
) -> dict[int, list[str]]:
    """
    只返回分组结果（不分配颜色），用于检查 / 调试分组。

    Returns
    -------
    dict[int, list[str]]
        group_id → [cell_type, ...]，每组首位是 dominant subtype。
    """
    groups, _, _, _ = _cluster_pipeline(coords, labels, n_major_groups)
    return groups


def color_h5ad(
    path: str,
    label_key: str = "cell_type",
    embedding_key: str = "X_umap",
    n_major_groups: int | None = None,
    palette: list[str] | None = None,
    return_groups: bool = False,
) -> dict[str, str] | tuple[dict[str, str], dict[int, list[str]]]:
    """
    从 h5ad 文件 (backed 模式) 直接生成配色，不加载表达矩阵。

    Parameters
    ----------
    path : str
        h5ad 文件路径。
    label_key : str
        obs 中的标签列名，默认 ``"cell_type"``。
    embedding_key : str
        obsm 中的降维 key，默认 ``"X_umap"``。
    n_major_groups, palette, return_groups : 同 ``assign_celltype_colors``。

    Returns
    -------
    dict[str, str]  或  (dict[str, str], dict[int, list[str]])
        label → hex color 映射。若 return_groups=True，额外返回分组。
    """
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    try:
        if embedding_key not in adata.obsm:
            avail = list(adata.obsm.keys())
            raise KeyError(
                f"Embedding '{embedding_key}' not found. Available: {avail}"
            )
        if label_key not in adata.obs.columns:
            avail = list(adata.obs.columns)
            raise KeyError(
                f"Label column '{label_key}' not found. Available: {avail}"
            )
        coords = np.asarray(adata.obsm[embedding_key][:, :2])
        labels = np.asarray(adata.obs[label_key])
    finally:
        adata.file.close()

    return assign_celltype_colors(
        coords, labels, n_major_groups, palette, return_groups=return_groups,
    )


# ══════════════════════════════════════════════════════════
#  可视化
# ══════════════════════════════════════════════════════════

def plot_umap(
    coords: np.ndarray,
    labels: np.ndarray,
    color_map: dict[str, str],
    groups: dict[int, list[str]] | None = None,
    *,
    # 布局
    figsize: tuple[float, float] = (18, 11),
    point_size: float = 1.2,
    # 质心标签
    show_group_labels: bool = True,
    label_fontsize: int = 10,
    label_fontweight: str = "bold",
    # 图例
    legend_fontsize: float = 7.5,
    legend_marker_size: float = 6,
    legend_columns: int | None = None,
    group_legend: bool = True,
    # 坐标轴
    xlabel: str = "UMAP 1",
    ylabel: str = "UMAP 2",
    # 输出
    save_path: str | None = None,
    dpi: int = 300,
    title: str | None = None,
):
    """
    绘制出版级 scatter plot。

    Features:
      - 分组图例（粗体组名 + 子项列表）
      - 大群质心标注
      - 稀有群渲染在最上层
      - 白色背景、无坐标轴刻度

    Parameters
    ----------
    coords, labels, color_map : 同 ``assign_celltype_colors``。
    groups : dict | None
        分组信息。``None`` 则内部自动推断。
    group_legend : bool
        True = 分组图例（带粗体标题），False = 平铺图例。
    xlabel, ylabel : str
        坐标轴标签，tSNE / PCA 时可改。
    其余参数控制样式和输出。

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe

    labels = np.asarray(labels)
    cell_counts = Counter(labels)

    if groups is None:
        groups = get_groups(coords, labels)

    sorted_gids = sorted(
        groups, key=lambda g: -sum(cell_counts[ct] for ct in groups[g])
    )
    group_names = {gid: _group_display_name(groups[gid]) for gid in sorted_gids}

    # ── scatter ──
    fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor="white")
    ax.set_facecolor("white")

    # 大群先画，稀有群后画（叠在上层）
    for ct in sorted(color_map, key=lambda c: cell_counts.get(c, 0), reverse=True):
        mask = labels == ct
        if not mask.any():
            continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=color_map[ct], s=point_size,
            edgecolors="none", rasterized=True, zorder=2,
        )

    # ── 质心标注 ──
    if show_group_labels:
        stroke = [pe.withStroke(linewidth=3, foreground="white")]
        for gid in sorted_gids:
            mask = np.isin(labels, groups[gid])
            if not mask.any():
                continue
            cx, cy = coords[mask].mean(axis=0)
            ax.text(
                cx, cy, group_names[gid],
                fontsize=label_fontsize, fontweight=label_fontweight,
                ha="center", va="center", color="#222222",
                zorder=5, path_effects=stroke,
            )

    # ── 图例 ──
    if group_legend:
        handles_labels: list[tuple] = []
        for gid in sorted_gids:
            handles_labels.append((
                Line2D([0], [0], marker="none", linestyle="none"),
                group_names[gid],
            ))
            for ct in groups[gid]:
                handles_labels.append((
                    Line2D(
                        [0], [0], marker="o", linestyle="none", color="w",
                        markerfacecolor=color_map[ct],
                        markersize=legend_marker_size, markeredgewidth=0,
                    ),
                    ct,
                ))
    else:
        handles_labels = [
            (
                Line2D(
                    [0], [0], marker="o", linestyle="none", color="w",
                    markerfacecolor=color_map[ct],
                    markersize=legend_marker_size, markeredgewidth=0,
                ),
                ct,
            )
            for gid in sorted_gids for ct in groups[gid]
        ]

    ncol = legend_columns or max(1, (len(handles_labels) + 24) // 25)
    leg = ax.legend(
        [h for h, _ in handles_labels],
        [l for _, l in handles_labels],
        fontsize=legend_fontsize,
        bbox_to_anchor=(1.01, 1), loc="upper left",
        frameon=False, ncol=ncol,
        handletextpad=0.3, labelspacing=0.35,
        columnspacing=1.0, borderpad=0,
    )
    # 组名加粗（不用 LaTeX，避免下划线等特殊字符问题）
    if group_legend:
        group_name_set = set(group_names.values())
        for text in leg.get_texts():
            if text.get_text() in group_name_set:
                text.set_fontweight("bold")

    # ── 样式 ──
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    if title:
        ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return fig, ax


def plot_palette(
    color_map: dict[str, str],
    groups: dict[int, list[str]],
    *,
    save_path: str | None = None,
):
    """
    按分组展示调色板预览。每行一个大群，★ 标记 dominant subtype。

    Returns
    -------
    (fig, ax)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    sorted_gids = sorted(groups, key=lambda g: -len(groups[g]))
    max_members = max(len(groups[g]) for g in sorted_gids)
    n_groups = len(sorted_gids)

    fig, ax = plt.subplots(
        1, 1,
        figsize=(min(20, max_members * 1.5 + 2), max(3, n_groups * 1.0)),
        facecolor="white",
    )
    ax.set_facecolor("white")

    for row, gid in enumerate(sorted_gids):
        members = groups[gid]
        y = -row * 1.1
        for col, ct in enumerate(members):
            color = color_map.get(ct, "#cccccc")
            x = col * 1.4
            rect = mpatches.FancyBboxPatch(
                (x, y), 1.2, 0.8, boxstyle="round,pad=0.06",
                facecolor=color, edgecolor="white", linewidth=1.5,
            )
            ax.add_patch(rect)
            display = f"★ {ct}" if col == 0 else ct
            text_color = "white" if _is_dark(color) else "#222222"
            ax.text(
                x + 0.6, y + 0.4, display,
                ha="center", va="center",
                fontsize=5.5, color=text_color, fontweight="bold",
            )
        ax.text(
            -0.3, y + 0.4, _group_display_name(members),
            ha="right", va="center",
            fontsize=9, fontweight="bold", color="#444444",
        )

    ax.set_xlim(-2, max_members * 1.4 + 0.5)
    ax.set_ylim(-n_groups * 1.1 + 0.2, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches="tight", facecolor="white")
    return fig, ax


# ══════════════════════════════════════════════════════════
#  内部工具
# ══════════════════════════════════════════════════════════

def _compute_centroids_fast(
    coords: np.ndarray, codes: np.ndarray, n_types: int,
) -> np.ndarray:
    """向量化质心计算：用整数编码 + bincount 代替逐类 Python 循环。"""
    counts = np.bincount(codes, minlength=n_types).astype(float)
    counts[counts == 0] = 1  # 避免除零
    sum_x = np.bincount(codes, weights=coords[:, 0], minlength=n_types)
    sum_y = np.bincount(codes, weights=coords[:, 1], minlength=n_types)
    return np.column_stack([sum_x / counts, sum_y / counts])


def _auto_determine_k(Z: np.ndarray, n_types: int) -> int:
    """在 dendrogram 合并距离中寻找自然切分点。

    策略：在搜索窗口内找到 relative gap 超过中位数 2 倍的第一个跳跃
    （从低 k 往高 k 扫描），避免总是被最后几次合并的巨大 gap 吸引。
    若无显著跳跃则回退到全局最大 gap。
    """
    if n_types <= 3:
        return n_types
    distances = Z[:, 2]
    gaps = np.diff(distances)
    rel_gaps = gaps / (distances[:-1] + 1e-10)

    k_max = max(15, n_types // 4)
    k_min = 3
    # 搜索窗口：对应 k 在 [k_min, k_max] 范围内的 gap 索引
    # gap index i 对应 k = n_types - i - 1
    idx_hi = len(rel_gaps) - k_min       # k=k_min 对应的 gap index
    idx_lo = max(0, len(rel_gaps) - k_max)  # k=k_max 对应的 gap index
    if idx_lo > idx_hi:
        idx_lo, idx_hi = idx_hi, idx_lo

    window = rel_gaps[idx_lo:idx_hi + 1]
    if len(window) == 0:
        return k_min

    median_gap = float(np.median(window))
    threshold = median_gap * 2.0

    # 从低 gap index（大 k）往高 gap index（小 k）扫描，
    # 找第一个超过阈值的跳跃 → 最细粒度的自然切分点
    best = idx_lo + int(np.argmax(window))  # 回退值：窗口内最大
    for i in range(idx_lo, idx_hi + 1):
        if rel_gaps[i] >= threshold:
            best = i
            break

    k = n_types - best - 1
    return max(k_min, min(k, k_max))


def _build_ordered_groups(
    unique_types: np.ndarray,
    major_labels: np.ndarray,
    type_counts: np.ndarray,
    Z: np.ndarray,
) -> dict[int, list[str]]:
    """构建分组映射，每组 dominant 排首位，其余按叶序。

    Parameters
    ----------
    type_counts : ndarray, shape (n_types,)
        每个 unique_type 的细胞计数（与 unique_types 对齐）。
    """
    # 建立 type → index 映射
    type_to_idx = {ct: i for i, ct in enumerate(unique_types)}

    groups: dict[int, list[str]] = {}
    for ct, ml in zip(unique_types, major_labels):
        groups.setdefault(int(ml), []).append(ct)

    leaf_order = leaves_list(Z)
    type_order = {unique_types[i]: r for r, i in enumerate(leaf_order)}

    for gid in groups:
        members = groups[gid]
        dominant = max(members, key=lambda ct: type_counts[type_to_idx[ct]])
        rest = sorted(
            [ct for ct in members if ct != dominant],
            key=lambda ct: type_order.get(ct, 0),
        )
        groups[gid] = [dominant] + rest
    return groups


def _common_prefix(strings: list[str]) -> str:
    """最长公共前缀，在词边界处截断。"""
    if not strings:
        return ""
    prefix = strings[0]
    for s in strings[1:]:
        while not s.startswith(prefix):
            prefix = prefix[:-1]
            if not prefix:
                return ""
    clean = prefix.rstrip(" _")
    if not clean:
        return ""
    if clean[-1].isdigit():
        for s in strings:
            if len(s) > len(clean) and s[len(clean)].isdigit():
                for i in range(len(clean) - 1, -1, -1):
                    if not clean[i].isdigit():
                        c = clean[:i + 1].rstrip(" _")
                        return c if len(c) >= 2 else ""
                return ""
    return clean


def _group_display_name(members: list[str]) -> str:
    """公共前缀 → 多数票首词 → dominant 名。"""
    prefix = _common_prefix(members)
    if len(prefix) >= 3:
        return prefix
    first_tokens = [re.split(r"[ _]", m)[0] for m in members if m]
    if first_tokens:
        most_common, count = Counter(first_tokens).most_common(1)[0]
        if count > len(members) * 0.4:
            return most_common
    return members[0]


def _is_dark(hex_color: str) -> bool:
    h = hex_color.lstrip("#")
    if len(h) == 3:
        h = h[0] * 2 + h[1] * 2 + h[2] * 2
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return 0.299 * r + 0.587 * g + 0.114 * b < 140
