# -------------------------
# Config / paths
# -------------------------
# addArchRGenome("hg38")

PATH_ARCHR_PROJ <- "./data/H3K27me3_coCnT_final3/"
PATH_MLIST_RDS <- "./data/list_matrix_imputed_gene_score.rds"
PATH_HDF5_STORE <- "./data/imputed_gene_score_h5"
PATH_HDF5_MANIFEST <- file.path(PATH_HDF5_STORE, "manifest.rds")
MARKER_CACHE_MAX <- 3
QUERY_CACHE_MAX <- 12
MAX_GENES_PER_QUERY <- 20

## Load data
# d_meta <- as.data.table(proj@cellColData)
# write_tsv(d_meta, "./data/cell_metadata.tsv")
d_meta <- fread("./data/cell_metadata.tsv")
umap_k27 <- fread("./data/table_umap_h3k27me3_final3.tsv")
setnames(umap_k27, c("UMAP1", "UMAP2"))
umap_k27[, cell_id := d_meta$cell_id]

message("[Shiny] Initializing marker matrix store...")

marker_store_mode <- "legacy_rds"
marker_manifest <- NULL
m_list <- NULL

if (file.exists(PATH_HDF5_MANIFEST)) {
    if (!requireNamespace("HDF5Array", quietly = TRUE)) {
        stop("HDF5 manifest found, but the 'HDF5Array' package is not installed.")
    }
    marker_manifest <- readRDS(PATH_HDF5_MANIFEST)
    if (!length(marker_manifest)) {
        stop("HDF5 manifest is empty.")
    }
    marker_store_mode <- "hdf5"
} else {
    message("[Shiny] HDF5 store not found. Falling back to monolithic RDS.")
    message("[Shiny] Loading m_list (imputed matrices)...")
    m_list <- readRDS(PATH_MLIST_RDS)
    stopifnot(length(m_list) > 0)
}

# -------------------------
# Marker colors
# -------------------------
marker_colors <- c(
    H3K27me3 = "#3B8FC4",
    H3K4me1 = "#FFD000",
    H3K4me2 = "#F39C12",
    H3K4me3 = "#E74C3C",
    H3K4me1_cooc = "#56B870",
    H3K4me2_cooc = "#D16BA5",
    H3K4me3_cooc = "#9B59B6"
)

default_markers <- c("H3K4me2_cooc", "H3K27me3", "H3K4me2")

# -------------------------
# Trajectory definitions
# -------------------------
trajectory_map <- list(
    "B cell" = list(name = "Bcell_Trajectory"),
    "Myeloid" = list(name = "Monocyte_Trajectory"),
    "Erythroid" = list(name = "Erythroid_Trajectory")
)

# -------------------------
# Load once (global)
# -------------------------
# message("[Shiny] Loading ArchR project...")
#
# library(ArchR)
# proj <- loadArchRProject(PATH_ARCHR_PROJ)
# umap_k27 <- as.data.table(getEmbedding(ArchRProj = proj, embedding = "UMAP_Harmony"))
# write_tsv(umap_k27, "./data/umap_coordinates.tsv")
#
# message("[Shiny] Ensuring trajectories exist...")
# for (nm in names(trajectory_map)) {
#     trj <- trajectory_map[[nm]]
#     proj <- addTrajectory(
#         ArchRProj = proj,
#         name = trj$name,
#         groupBy = "Clusters_plus",
#         trajectory = trj$clusters,
#         embedding = "UMAP_Harmony",
#         force = TRUE
#     )
# }

marker_cache <- make_lru_cache(MARKER_CACHE_MAX)

if (identical(marker_store_mode, "hdf5")) {
    available_markers <- intersect(names(marker_manifest), names(marker_colors))
} else {
    available_markers <- intersect(names(m_list), names(marker_colors))
}

if (length(available_markers) == 0) {
    stop("No markers available after initializing the matrix store.")
}

available_genes <- rownames(get_marker_matrix(available_markers[[1]]))
if (is.null(available_genes)) {
    stop("Marker matrices must have rownames as genes.")
}

gene_index <- stats::setNames(seq_along(available_genes), available_genes)
default_markers <- intersect(default_markers, available_markers)
if (length(default_markers) == 0) {
    default_markers <- available_markers[1]
}
