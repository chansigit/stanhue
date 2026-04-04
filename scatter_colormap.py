"""
umap_colormap.py — 通用单细胞 UMAP 层级化自动配色工具
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
  >>> from umap_colormap import assign_celltype_colors, plot_umap
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

def assign_celltype_colors(
    coords: np.ndarray,
    labels: np.ndarray,
    n_major_groups: int | None = None,
    palette: list[str] | None = None,
) -> dict[str, str]:
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

    Returns
    -------
    dict[str, str]
        cell_type → hex color 映射。
    """
    coords, labels = _validate_inputs(coords, labels)

    if palette is None:
        palette = PAIRED_PALETTE
    n_pal = len(palette)

    unique_types = np.unique(labels)
    n_types = len(unique_types)

    # 边界情况
    if n_types == 0:
        return {}
    if n_types == 1:
        return {unique_types[0]: palette[0]}
    if n_types <= n_pal and n_major_groups is None:
        return {ct: palette[i] for i, ct in enumerate(unique_types)}

    # Step 1–3: 聚类 + 分组 + 组内排序
    centroids = _compute_centroids(coords, labels, unique_types)
    Z = linkage(centroids, method="ward")

    if n_major_groups is None:
        n_major_groups = _auto_determine_k(Z, n_types)
    n_major_groups = min(n_major_groups, n_types)

    major_labels = fcluster(Z, t=n_major_groups, criterion="maxclust")
    groups = _build_ordered_groups(unique_types, major_labels, labels, Z)

    # Step 4: 按总细胞数排序，分配偏移
    cell_counts = Counter(labels)
    sorted_gids = sorted(
        groups, key=lambda g: -sum(cell_counts[ct] for ct in groups[g])
    )
    n_groups = len(sorted_gids)

    if n_groups <= n_pal // 2:
        offsets = [i * 2 for i in range(n_groups)]           # 步长 2
    else:
        offsets = list(range(0, n_pal, 2)) + list(range(1, n_pal, 2))
        offsets = offsets[:n_groups]

    # Step 5: 组内顺延取色
    color_map: dict[str, str] = {}
    for gid, offset in zip(sorted_gids, offsets):
        for i, ct in enumerate(groups[gid]):
            color_map[ct] = palette[(offset + i) % n_pal]

    return color_map


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
    labels = np.asarray(labels)
    unique_types = np.unique(labels)
    centroids = _compute_centroids(coords, labels, unique_types)
    Z = linkage(centroids, method="ward")
    if n_major_groups is None:
        n_major_groups = _auto_determine_k(Z, len(unique_types))
    n_major_groups = min(n_major_groups, len(unique_types))
    major_labels = fcluster(Z, t=n_major_groups, criterion="maxclust")
    return _build_ordered_groups(unique_types, major_labels, labels, Z)


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
                f"$\\bf{{{group_names[gid]}}}$",
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
    ax.legend(
        [h for h, _ in handles_labels],
        [l for _, l in handles_labels],
        fontsize=legend_fontsize,
        bbox_to_anchor=(1.01, 1), loc="upper left",
        frameon=False, ncol=ncol,
        handletextpad=0.3, labelspacing=0.35,
        columnspacing=1.0, borderpad=0,
    )

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

def _compute_centroids(
    coords: np.ndarray, labels: np.ndarray, unique_types: np.ndarray,
) -> np.ndarray:
    centroids = np.zeros((len(unique_types), 2))
    for i, ct in enumerate(unique_types):
        centroids[i] = coords[labels == ct].mean(axis=0)
    return centroids


def _auto_determine_k(Z: np.ndarray, n_types: int) -> int:
    """在 dendrogram 合并距离中找最大相对间距跳跃。"""
    if n_types <= 3:
        return n_types
    distances = Z[:, 2]
    gaps = np.diff(distances)
    rel_gaps = gaps / (distances[:-1] + 1e-10)
    start = max(0, len(rel_gaps) - min(25, n_types - 1))
    best = start + np.argmax(rel_gaps[start:])
    return max(3, min(n_types - best - 1, 15))


def _build_ordered_groups(
    unique_types: np.ndarray,
    major_labels: np.ndarray,
    all_labels: np.ndarray,
    Z: np.ndarray,
) -> dict[int, list[str]]:
    """构建分组映射，每组 dominant 排首位，其余按叶序。"""
    groups: dict[int, list[str]] = {}
    for ct, ml in zip(unique_types, major_labels):
        groups.setdefault(int(ml), []).append(ct)

    cell_counts = Counter(all_labels)
    leaf_order = leaves_list(Z)
    type_order = {unique_types[i]: r for r, i in enumerate(leaf_order)}

    for gid in groups:
        members = groups[gid]
        dominant = max(members, key=lambda ct: cell_counts[ct])
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
    r, g, b = (int(hex_color[i:i + 2], 16) for i in (1, 3, 5))
    return 0.299 * r + 0.587 * g + 0.114 * b < 140
