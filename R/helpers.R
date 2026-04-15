make_lru_cache <- function(max_size) {
    cache_env <- new.env(parent = emptyenv())
    cache_order <- character(0)

    list(
        exists = function(key) {
            exists(key, envir = cache_env, inherits = FALSE)
        },
        get = function(key) {
            cache_order <<- c(setdiff(cache_order, key), key)
            get(key, envir = cache_env, inherits = FALSE)
        },
        set = function(key, value) {
            cache_order <<- c(setdiff(cache_order, key), key)
            assign(key, value, envir = cache_env)

            if (length(cache_order) > max_size) {
                evict <- cache_order[[1]]
                if (exists(evict, envir = cache_env, inherits = FALSE)) {
                    rm(list = evict, envir = cache_env)
                }
                cache_order <<- cache_order[-1]
            }
        }
    )
}

get_marker_matrix <- function(marker) {
    if (marker_cache$exists(marker)) {
        return(marker_cache$get(marker))
    }

    mat <- if (identical(marker_store_mode, "hdf5")) {
        spec <- marker_manifest[[marker]]
        if (is.null(spec)) {
            stop(sprintf("Marker '%s' is not available in the HDF5 store.", marker))
        }
        HDF5Array::HDF5Array(filepath = spec$filepath, name = spec$name)
    } else {
        m_list[[marker]]
    }

    if (is.null(mat)) {
        stop(sprintf("Marker '%s' could not be loaded.", marker))
    }

    marker_cache$set(marker, mat)
    mat
}

normalize_gene_matrix <- function(mat, method = c("z", "minmax")) {
    method <- match.arg(method)

    if (method == "z") {
        mu <- rowMeans(mat, na.rm = TRUE)
        sdv <- apply(mat, 1, sd, na.rm = TRUE)
        sdv[sdv == 0 | is.na(sdv)] <- 1
        mat <- sweep(mat, 1, mu, "-")
        mat <- sweep(mat, 1, sdv, "/")
    }

    if (method == "minmax") {
        rmin <- apply(mat, 1, min, na.rm = TRUE)
        rmax <- apply(mat, 1, max, na.rm = TRUE)
        denom <- rmax - rmin
        denom[denom == 0 | is.na(denom)] <- 1
        mat <- sweep(mat, 1, rmin, "-")
        mat <- sweep(mat, 1, denom, "/")
    }

    mat
}

parse_genes_from_text <- function(txt) {
    g <- unlist(strsplit(txt, "[,\n\r\t ]+"))
    g <- toupper(trimws(g))
    g <- g[nzchar(g)]
    unique(g)
}

sanitize_filename_component <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", x)
}

make_export_filename <- function(prefix, payload, ext) {
    sprintf(
        "%s_%s_%s.%s",
        prefix,
        sanitize_filename_component(payload$lineage),
        sanitize_filename_component(payload$geneset_label),
        ext
    )
}

make_status_message <- function(prefix, geneset_label, lineage, marker_count, point_count) {
    sprintf(
        "%s: %s | lineage: %s | markers: %d | points: %d",
        prefix,
        geneset_label,
        lineage,
        marker_count,
        point_count
    )
}

scientific_theme <- function(base_size = 14) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            plot.title = element_text(face = "bold", size = rel(1.15), margin = margin(b = 8)),
            plot.subtitle = element_text(color = "#52606D", size = rel(0.92), margin = margin(b = 12)),
            plot.caption = element_text(color = "#6B7280", size = rel(0.82), hjust = 1),
            axis.title = element_text(face = "bold", color = "#102A43"),
            axis.text = element_text(color = "#243B53"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "#D9E2EC", linewidth = 0.35),
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            legend.text = element_text(color = "#243B53"),
            strip.text = element_text(face = "bold", color = "#102A43"),
            strip.background = element_rect(fill = "#F0F4F8", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
}

plot_subtitle <- function(lineage_label, geneset_label, marker_count) {
    sprintf(
        "%s lineage | %s | %d marker%s",
        lineage_label,
        geneset_label,
        marker_count,
        if (marker_count == 1) "" else "s"
    )
}

make_geneset_score_dt <- function(
    genes,
    markers,
    allow_missing = TRUE,
    normalize = FALSE,
    norm_method = c("z", "minmax")
) {
    norm_method <- match.arg(norm_method)
    genes <- toupper(trimws(genes))
    genes <- genes[nzchar(genes)]
    genes <- unique(genes)
    if (length(genes) == 0) {
        stop("No genes provided.")
    }

    markers <- intersect(markers, available_markers)
    if (length(markers) == 0) {
        stop("No valid markers selected.")
    }

    present <- intersect(genes, available_genes)
    missing <- setdiff(genes, available_genes)

    if (!allow_missing && length(missing) > 0) {
        stop(sprintf(
            "Missing genes (not in matrix rownames): %s",
            paste(missing, collapse = ", ")
        ))
    }
    if (length(present) == 0) {
        stop("None of the input genes are found in available_genes.")
    }

    row_idx <- unname(gene_index[present])
    row_idx <- row_idx[!is.na(row_idx)]

    dt_list <- lapply(markers, function(mk) {
        m <- get_marker_matrix(mk)
        if (length(row_idx) == 0) {
            return(NULL)
        }

        sub <- as.matrix(m[row_idx, , drop = FALSE])

        if (normalize) {
            sub <- normalize_gene_matrix(sub, method = norm_method)
        }

        v <- colMeans(sub, na.rm = TRUE)

        data.table(
            cell_id = names(v),
            gene_score = as.numeric(v),
            marker = mk
        )
    })

    d_out <- rbindlist(dt_list, fill = TRUE)
    if (nrow(d_out) == 0) {
        stop("No data produced for the selected markers/genes.")
    }

    list(
        dt = d_out,
        present = present,
        missing = missing
    )
}

make_gene_score_dt <- function(gene, markers) {
    gene <- toupper(trimws(gene))
    if (!nzchar(gene)) {
        stop("Gene is empty.")
    }

    if (!(gene %in% available_genes)) {
        stop(sprintf("Gene '%s' not found in m_list rownames.", gene))
    }

    markers <- intersect(markers, available_markers)
    if (length(markers) == 0) {
        stop("No valid markers selected.")
    }

    gene_idx <- unname(gene_index[gene])
    if (is.na(gene_idx)) {
        stop(sprintf("Gene '%s' not found in m_list rownames.", gene))
    }

    lapply(markers, function(mk) {
        m <- get_marker_matrix(mk)
        v <- m[gene_idx, ]
        data.table(
            cell_id = names(v),
            gene_score = as.numeric(v),
            marker = mk
        )
    }) %>% rbindlist(fill = TRUE)
}

make_pseudotime_dt <- function(lineage_label) {
    trj_col <- trajectory_map[[lineage_label]]$name
    if (!(trj_col %in% colnames(d_meta))) {
        stop(sprintf("Trajectory column '%s' not found in ArchR cellColData.", trj_col))
    }

    d_meta[, .(cell_id, pseudotime = get(trj_col))]
}

plot_gene_vs_pseudotime <- function(d_plot, gene, lineage_label, selected_markers) {
    d_plot[, marker := factor(marker, levels = selected_markers)]

    ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = marker)) +
        geom_point(alpha = 0.38, size = 0.9, stroke = 0) +
        geom_smooth(method = "gam", se = FALSE, linewidth = 1.05, span = 0.08, n = 300) +
        scale_color_manual(values = marker_colors, drop = FALSE) +
        labs(
            title = sprintf("%s Across Pseudotime", gene),
            subtitle = plot_subtitle(lineage_label, gene, length(selected_markers)),
            x = "Pseudotime",
            y = "Gene Set Score",
            color = "Marker"
        ) +
        scientific_theme()
}
