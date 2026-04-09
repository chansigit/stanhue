---
name: stanhue
description: >
  Use ONLY when the user explicitly mentions "stanhue" by name.
  Do NOT auto-trigger on generic color/palette/UMAP keywords.
  Do NOT use for gene harmonization (that is stangene) or format
  conversion (that is stanobj).
version: 1.0.0
allowed-tools: [Bash, Read, Write, Glob, Grep]
---

# Skill: Hierarchical Auto-Coloring for Scatter Plots

Generate a `{label: hex_color}` palette for any scatter plot with categorical
labels, based on spatial relationships in 2D coordinates. This skill **only
produces the color mapping** — it does NOT plot or visualize anything.
The user is responsible for applying the returned palette to their own
plotting code.

Typical use cases: single-cell UMAP/tSNE with cell type labels, but equally
applicable to any 2D scatter plot where categories need coherent coloring.

Uses Ward hierarchical clustering on category centroids to infer groups,
then distributes colors from a repeating Paired palette with designed
offsets so that:
- **Distant lineages** get different hue families (blue vs red vs green)
- **Related subtypes** get adjacent shades within the same hue family
- **Dominant subtypes** (most cells) anchor each group's representative color

## Algorithm Overview

```
UMAP coords + cell type labels
        │
        ▼
  1. Compute centroid per cell type
        │
        ▼
  2. Ward hierarchical clustering on centroids
     → auto-determine k (relative gap method, k ∈ [3,15])
        │
        ▼
  3. Build ordered groups:
     - dominant subtype (most cells) → position 0
     - rest sorted by dendrogram leaf order
        │
        ▼
  4. Assign palette offsets to groups:
     - groups sorted by total cell count (descending)
     - offset step = 2 in Paired palette (0→2→4→6→8→10)
     - if >6 groups, fill odd positions (1→3→5→7→9→11)
        │
        ▼
  5. Each group walks palette from its offset (mod palette_len)
        │
        ▼
  Output: { cell_type: "#hex_color" }
```

## Prerequisites

The scripts are bundled with this skill. Find the script directory:

```bash
# Plugin mode (preferred)
STANHUE_SCRIPTS=$(find ~/.claude/plugins/cache -path "*/stanhue/skills/stanhue/scripts" -type d 2>/dev/null | head -1)
# Fallback: local skill
[ -z "$STANHUE_SCRIPTS" ] && STANHUE_SCRIPTS="$HOME/.claude/skills/stanhue/scripts"
ls "$STANHUE_SCRIPTS/scatter_colormap.py" "$STANHUE_SCRIPTS/scatter_colormap.R"
```

**Python dependencies:** `numpy`, `scipy`
**R dependencies:** `stats` (base R only)

## Input Validation

Before running, validate the user's inputs:

1. **coords** must be a numeric matrix/array with shape `(n_cells, 2)`.
   - Reject if not 2 columns, contains NaN/Inf, or is empty.
2. **labels** must be a 1D array/vector of strings with length `n_cells`.
   - Reject if length mismatches coords rows.
   - Reject if all labels are identical (nothing to color-differentiate).
3. **At least 2 unique cell types** are needed for the clustering to be meaningful.
   - 1 unique type → return single color, no clustering needed.
   - ≤12 unique types and no manual k → simple sequential palette assignment.

```python
# Python validation example
import numpy as np

def validate_inputs(coords, labels):
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    errors = []
    if coords.ndim != 2 or coords.shape[1] != 2:
        errors.append(f"coords must be (n, 2), got shape {coords.shape}")
    if np.any(~np.isfinite(coords)):
        errors.append("coords contains NaN or Inf values")
    if len(labels) != coords.shape[0]:
        errors.append(f"labels length ({len(labels)}) != coords rows ({coords.shape[0]})")
    if len(np.unique(labels)) == 0:
        errors.append("no cell type labels provided")
    if errors:
        raise ValueError("Input validation failed:\n" + "\n".join(f"  - {e}" for e in errors))
    return coords, labels
```

```r
# R validation example
validate_inputs <- function(coords, labels) {
  coords <- as.matrix(coords)
  labels <- as.character(labels)
  errors <- character(0)
  if (ncol(coords) != 2)
    errors <- c(errors, sprintf("coords must have 2 columns, got %d", ncol(coords)))
  if (any(!is.finite(coords)))
    errors <- c(errors, "coords contains NaN or Inf values")
  if (length(labels) != nrow(coords))
    errors <- c(errors, sprintf("labels length (%d) != coords rows (%d)",
                                length(labels), nrow(coords)))
  if (length(unique(labels)) == 0)
    errors <- c(errors, "no cell type labels provided")
  if (length(errors) > 0)
    stop("Input validation failed:\n", paste("  -", errors, collapse = "\n"))
  list(coords = coords, labels = labels)
}
```

## Usage — Python

```python
import subprocess, sys
scripts_dir = subprocess.check_output(
    'find ~/.claude/plugins/cache -path "*/stanhue/skills/stanhue/scripts" -type d 2>/dev/null | head -1',
    shell=True, text=True).strip()
if not scripts_dir:
    scripts_dir = os.path.expanduser("~/.claude/skills/stanhue/scripts")
sys.path.insert(0, scripts_dir)
from scatter_colormap import assign_celltype_colors, get_groups, color_h5ad

# ---- Fastest: directly from h5ad (backed mode, skips expression matrix) ----
color_map = color_h5ad("data.h5ad", label_key="cell_type", embedding_key="X_umap")

# ---- From arrays ----
color_map = assign_celltype_colors(
    coords=adata.obsm["X_umap"],   # (n_cells, 2) array
    labels=adata.obs["cell_type"], # (n_cells,) labels
)
# color_map: {"CD4 Naive": "#a6cee3", "CD8 TEM": "#1f78b4", ...}

# Optional: inspect grouping
groups = get_groups(adata.obsm["X_umap"], adata.obs["cell_type"])
```

## Usage — R

```r
# Use the STANHUE_SCRIPTS path found in Prerequisites
source(file.path(Sys.getenv("STANHUE_SCRIPTS", "~/.claude/skills/stanhue/scripts"), "scatter_colormap.R"))

# Core: get color mapping
color_map <- assign_celltype_colors(
  coords = Embeddings(seurat_obj, "umap"),  # matrix (n x 2)
  labels = seurat_obj$cell_type             # character vector
)
# Named vector: c("CD4 Naive" = "#a6cee3", ...)

# Convenience wrappers for Seurat / SCE objects:
color_map <- color_seurat(seurat_obj, reduction = "umap", group_by = "cell_type")
color_map <- color_sce(sce_obj, dimred = "UMAP", col_name = "cell_type")

# Inspect grouping
groups <- get_groups(coords, labels)
```

## Custom Palette

The default `PAIRED_PALETTE` (ColorBrewer Paired, 12 colors) works well for most
cases. To use a custom palette, pass any ordered list of hex colors:

```python
my_palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"]
color_map = assign_celltype_colors(coords, labels, palette=my_palette)
```

The offset logic adapts automatically: step size = `len(palette) // n_groups`
(minimum 1), so groups always start on distinct colors when possible.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `coords` | required | 2D embedding coordinates, shape (n, 2) |
| `labels` | required | Cell type labels, length n |
| `n_major_groups` | `None`/auto | Number of lineage groups (auto uses relative gap) |
| `palette` | `PAIRED_PALETTE` | Ordered hex color list |

## When Results Need Manual Adjustment

- If two visually distinct clusters share a color, increase `n_major_groups`.
- If related subtypes get unrelated colors, decrease `n_major_groups`.
- For >30 cell types, consider a larger custom palette (e.g., 20-color set).

## Scope

- This skill **only generates the palette** (`{label: hex_color}` mapping).
- It does **NOT** plot, render, or visualize anything.
- The user applies the returned palette to their own plotting code.

## Important

- The algorithm is **deterministic** — same input always produces same colors.
- Works with any 2D coordinates + categorical labels — not limited to single-cell
  data. Any scatter plot with category labels can use this algorithm.
- Does NOT require lineage or domain annotations — infers structure from spatial layout.
- For Seurat objects in R, use the `color_seurat()` convenience wrapper.
- Report the grouping structure to the user so they can verify it makes sense.
