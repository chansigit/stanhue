"""stanhue — hierarchical, overlap-aware auto-coloring for scatter plots with categorical labels.

>>> from stanhue import assign_celltype_colors
>>> colors = assign_celltype_colors(coords, labels)
"""

from .scatter_colormap import (
    PAIRED_PALETTE,
    assign_celltype_colors,
    color_h5ad,
    get_groups,
    plot_palette,
    plot_umap,
)

__version__ = "1.1.1"

__all__ = [
    "PAIRED_PALETTE",
    "assign_celltype_colors",
    "color_h5ad",
    "get_groups",
    "plot_palette",
    "plot_umap",
    "__version__",
]
