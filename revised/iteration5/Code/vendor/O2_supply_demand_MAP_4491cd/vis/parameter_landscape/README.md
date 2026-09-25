# Parameter-landscape visualization

| File | Responsibility |
|---|---|
| `visualize_parameter_landscape.R` | Reads analysis tables and renders PCA/t-SNE/UMAP, cluster, and contribution figures. |
| `parameter_landscape_visualization_utils.R` | Plot constructors, palettes, and figure-save helpers; no model or fit-object reads. |

All plotted values must already exist in the analysis tables.
