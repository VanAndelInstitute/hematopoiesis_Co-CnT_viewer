# Hematopoiesis Single cell CoCut&Tag histone modification viewer

This repository contains a Shiny app for exploring single-cell CoCut&Tag bone marrow data across lineage pseudotime trajectories. The app projects single cell CoCut&Tag gene score onto pseudotime and UMAP space, and supports export of plots and raw tables.

<img width="800" height="484" alt="demo" src="https://github.com/user-attachments/assets/fe5df9d0-58c1-4156-a289-32ff4a2981aa" />

## Features

- Query one gene or a small gene set.
- Compare multiple histone marks in the same view.
- Switch between lineage trajectories: `B cell`, `Myeloid`, and `Erythroid`.
- Render synchronized pseudotime and UMAP plots.
- Export figures as PDF and merged data as TSV.

## Repository Layout

- `app.R`: main Shiny entry point.
- `dev_app.R`: development launcher with autoreload enabled.
- `R/`
  - `globals.R`: paths, data loading, marker colors, lineage definitions.
  - `helpers.R`: score calculation, caching, plotting helpers.
  - `ui.R`: user interface and styling.
  - `server.R`: reactive logic and download handlers.
- `data/`: app assets and input data (TODO: consider put large data on Zenodo and load on demand).

## Data Inputs

The app expects these local files and directories:

- `data/cell_metadata.tsv`
- `data/table_umap_h3k27me3_final3.tsv`
- `data/imputed_gene_score_h5/manifest.rds`
- `data/imputed_gene_score_h5/*.h5`

Static asset:

- `data/fig_dataset_diagram.png` is served directly by the app for the dataset panel.

## R Dependencies

Required packages:

- `shiny`
- `data.table`
- `ggplot2`
- `magrittr`
- `HDF5Array`

Example install:

```r
install.packages(c("shiny", "data.table", "ggplot2", "magrittr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("HDF5Array")
```

## Run The App

From the repository root:

```bash
Rscript app.R
```

For development with autoreload:

```bash
Rscript dev_app.R
```

