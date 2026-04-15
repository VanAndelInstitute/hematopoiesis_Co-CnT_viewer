#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(shiny)
    library(data.table)
    library(ggplot2)
    library(magrittr)
})

# -------------------------
# Config / paths
# -------------------------
# addArchRGenome("hg38")

PATH_ARCHR_PROJ <- "./data/H3K27me3_coCnT_final3/"
PATH_MLIST_RDS  <- "./data/list_matrix_imputed_gene_score.rds"

## Load data
# d_meta <- as.data.table(proj@cellColData)
# write_tsv(d_meta, "./data/cell_metadata.tsv")
d_meta <- fread("./data/cell_metadata.tsv")
dim(d_meta)
umap_k27 = fread("./data/table_umap_h3k27me3_final3.tsv")

message("[Shiny] Loading m_list (imputed matrices)...")
m_list <- readRDS(PATH_MLIST_RDS)
stopifnot(length(m_list) > 0)



# -------------------------
# Marker colors (yours)
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
    "B cell"    = list(name = "Bcell_Trajectory"),
    "Monocyte"  = list(name = "Monocyte_Trajectory"),
    "Erythroid" = list(name = "Erythroid_Trajectory")
)

# -------------------------
# Load once (global)
# -------------------------
# message("[Shiny] Loading ArchR project...")

# library(ArchR)
# proj <- loadArchRProject(PATH_ARCHR_PROJ)
# umap_k27 = as.data.table( getEmbedding( ArchRProj = proj, embedding = "UMAP_Harmony"))
# write_tsv(umap_k27, "./data/umap_coordinates.tsv")


# message("[Shiny] Ensuring trajectories exist...")
# for (nm in names(trajectory_map)) {
#     trj <- trajectory_map[[nm]]
#     proj <- addTrajectory(
# 	ArchRProj  = proj,
# 	name       = trj$name,
# 	groupBy    = "Clusters_plus",
# 	trajectory = trj$clusters,
# 	embedding  = "UMAP_Harmony",
# 	force      = TRUE
#     )
# }

available_markers <- intersect(names(m_list), names(marker_colors))
available_genes <- rownames(m_list[[1]])
if (is.null(available_genes)) stop("m_list[[1]] must have rownames as genes.")

default_markers <- intersect(default_markers, available_markers)
if (length(default_markers) == 0) default_markers <- available_markers[1]

# -------------------------
# Helpers
# -------------------------
normalize_gene_matrix <- function(mat, method = c("z", "minmax")) {
    method <- match.arg(method)

    if (method == "z") {
	# z-score per gene (row)
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
    # one gene per line, but also tolerate commas/spaces
    g <- unlist(strsplit(txt, "[,\n\r\t ]+"))
    g <- toupper(trimws(g))
    g <- g[nzchar(g)]
    unique(g)
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
    if (length(genes) == 0) stop("No genes provided.")

    markers <- intersect(markers, names(m_list))
    if (length(markers) == 0) stop("No valid markers selected.")

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

    # For each marker: avg across genes per cell (colMeans over selected rows)
    dt_list <- lapply(markers, function(mk) {
	m <- m_list[[mk]]
	idx <- intersect(present, rownames(m))
	if (length(idx) == 0) return(NULL)


	sub <- as.matrix(m[idx, , drop = FALSE])

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
    if (nrow(d_out) == 0) stop("No data produced for the selected markers/genes.")

    list(
	dt = d_out,
	present = present,
	missing = missing
    )
}

make_gene_score_dt <- function(gene, markers) {
    gene <- toupper(trimws(gene))
    if (!nzchar(gene)) stop("Gene is empty.")

    if (!(gene %in% available_genes)) {
	stop(sprintf("Gene '%s' not found in m_list rownames.", gene))
    }

    markers <- intersect(markers, names(m_list))
    if (length(markers) == 0) stop("No valid markers selected.")

    lapply(markers, function(mk) {
	m <- m_list[[mk]]
	if (!(gene %in% rownames(m))) return(NULL)
	v <- m[gene, ]
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
    d_meta[, .(cell_id, pseudotime = get(trj_col))] #%>% na.omit()
}

plot_gene_vs_pseudotime <- function(d_plot, gene, lineage_label, selected_markers) {
    d_plot[, marker := factor(marker, levels = selected_markers)]

    ggplot(d_plot, aes(x = pseudotime, y = gene_score, color = marker)) +
	geom_point(alpha = 0.5, size = 1) +
	geom_smooth(method = "gam", se = FALSE, span = 0.05, n = 200) +
	scale_color_manual(values = marker_colors, drop = FALSE) +
	theme_classic() +
	theme(
	    legend.position = "right",
	    axis.title = element_blank()
	    ) +
	ggtitle(sprintf("%s — %s pseudotime", gene, lineage_label))
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
    titlePanel("Gene-set (avg) score vs lineage pseudotime"),

    sidebarLayout(
	sidebarPanel(
	    selectInput(
		"lineage", "Lineage / trajectory",
		choices = names(trajectory_map),
		selected = "B cell"
		),

	    # textInput(
	    # "gene", "Gene",
	    # value = "DNMT3B",
	    # placeholder = "Type gene symbol, e.g. DNMT3B"
	    # ),

	    textAreaInput(
		"genes", "Genes (one per line)",
		value = "DNMT3B\nDNMT3A\nTET1",
		placeholder = "DNMT3B\nDNMT3A\nTET1",
		rows = 6,
		width = "100%"
		),

	    checkboxInput(
		"normalize_genes",
		"Normalize each gene before averaging",
		value = FALSE
		),

	    conditionalPanel(
		condition = "input.normalize_genes == true",
		selectInput(
		    "norm_method",
		    "Normalization method",
		    choices = c(
			"zscore (per gene)" = "z",
			"min-max (0–1)"     = "minmax"
			),
		    selected = "z"
		)
		),

	    checkboxInput("allow_missing_genes", "Allow missing genes (ignore)", value = TRUE),


	    selectizeInput(
		"markers", "Markers",
		choices = available_markers,
		selected = default_markers,
		multiple = TRUE,
		options = list(plugins = list("remove_button"), maxOptions = 5000)
		),

	    checkboxInput("facet_markers", "Facet by marker", value = FALSE),

	    actionButton("submit", "Submit / Plot", class = "btn-primary"),
	    hr(),
	    verbatimTextOutput("status")
	    ),

	mainPanel(
	    fluidRow(
		column(6, downloadButton("download_pseudotime_pdf", "Download pseudotime PDF")),
		column(6, downloadButton("download_umap_pdf", "Download UMAP PDF")),
		column(4, downloadButton("download_raw_table", "Download raw table (TSV)"))

		),
	    plotOutput("p_pseudotime", height = "650px"),
	    plotOutput("p_umap", height = "650px")
	)
    )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {

    status_msg <- reactiveVal("Choose gene/lineage/markers, then click Submit / Plot.")

    # eventReactive: only recompute when submit is clicked
    plot_payload <- eventReactive(input$submit, {
	lineage <- isolate(input$lineage)

	normalize <- isolate(input$normalize_genes)
	norm_method <- isolate(input$norm_method)

	print("OK6---------------------------------------")

	genes_txt <- isolate(input$genes)
	allow_missing_genes <- isolate(input$allow_missing_genes)
	print("OK5---------------------------------------")

	markers <- isolate(input$markers)
	facet <- isolate(input$facet_markers)

	print("OK4---------------------------------------")
	# print(lineage)
	# print(gene)
	# print(markers)
	# print(facet)

	# lineage = "B cell"
	# gene = "PAX5"
	# markers = c("H3K4me2_cooc", "H3K27me3")
	# facet = FALSE

	# genes_txt = "PAX5\nRAG1\nRAG2"
	# markers = c("H3K4me2_cooc", "H3K27me3")
	# lineage = "B cell"
	# allow_missing_genes = TRUE



	names(umap_k27) = c("UMAP1", "UMAP2")
	umap_k27[, cell_id := d_meta$cell_id]
	print("OK3---------------------------------------")

	genes = parse_genes_from_text(genes_txt)
	selected_markers <- intersect(markers, available_markers)
	print("OK2---------------------------------------")

	if (length(genes) == 0) stop("Please input at least one gene (one per line).")
	if (length(selected_markers) == 0) stop("Please select at least one marker.")

	print("OK1---------------------------------------")

	pt <- make_pseudotime_dt(lineage)

	# gs <- make_gene_score_dt(gene, selected_markers)
	gs_obj = make_geneset_score_dt(
	    genes, 
	    selected_markers, 
	    allow_missing = allow_missing_genes,
	    normalize = normalize,
	    norm_method = norm_method
	)

	gs = gs_obj$dt

	d <- merge(gs, pt, by = "cell_id")
	d <- merge(d, umap_k27, by = "cell_id", all.y = T)
	# d <- d[is.finite(pseudotime) & is.finite(gene_score)]
	if (nrow(d) == 0) stop("No data after merging gene scores with pseudotime.")


	geneset_label <- sprintf("Avg(%d genes)", length(gs_obj$present))
	if (length(gs_obj$missing) > 0) {
	    geneset_label <- sprintf(
		"%s | missing:%d",
		geneset_label, length(gs_obj$missing)
	    )
	}


	if (length(gs_obj$missing) > 0) {
	    status_msg(paste0(
		    status_msg(),
		    "\nUsed: ", paste(gs_obj$present, collapse = ", "),
		    "\nMissing: ", paste(gs_obj$missing, collapse = ", ")
		    ))
	}

	status_msg(sprintf(
		"OK: %s | lineage: %s | markers: %d | points: %d",
		geneset_label, lineage, length(unique(d$marker)), nrow(d)
		))

	list(
	    d = d,
	    genes = gs_obj$present,
	    missing = gs_obj$missing,
	    geneset_label = geneset_label,
	    lineage = lineage,
	    selected_markers = selected_markers,
	    facet = facet
	)
    }, ignoreInit = TRUE)


    # ---- build ggplot objects once (reused for renderPlot + download) ----
    pseudotime_plot <- reactive({
	payload <- plot_payload()
	req(payload)
	d <- payload$d

	if (isTRUE(payload$facet)) {
	    ggplot(d, aes(x = pseudotime, y = gene_score)) +
		geom_point(alpha = 0.5, size = 1) +
		geom_smooth(method = "gam", se = FALSE, span = 0.05, n = 200) +
		facet_wrap(~marker, scales = "free_y") +
		theme_classic() +
		theme(axis.title = element_blank()) +
		ggtitle(sprintf("%s — %s pseudotime", payload$geneset_label, payload$lineage))
	} else {
	    plot_gene_vs_pseudotime(
		d, payload$geneset_label, payload$lineage, payload$selected_markers
	    )
	}
    })

    umap_plot <- reactive({
	payload <- plot_payload()
	req(payload)

	d <- payload$d[, .(UMAP1, UMAP2, gene_score, marker)] %>% na.omit()

	high_color = if (length(payload$selected_markers) == 1) {
	    marker_colors[payload$selected_markers]
	} else {
	    "red"
	}

	ggplot(d, aes(x = -UMAP1, y = UMAP2, color = gene_score)) +
	    geom_point(alpha = 0.5, size = 1) +
	    facet_wrap(~marker, scale = "free") +
	    scale_color_continuous(low = "lightgrey", high = high_color) +
	    theme_classic() +
	    theme(axis.title = element_blank()) +
	    ggtitle(sprintf("%s score on UMAP — %s lineage", payload$geneset_label, payload$lineage))
    })

    # ---- renderPlot uses the same ggplot objects ----
    output$p_pseudotime <- renderPlot({
	p <- pseudotime_plot()
	print(p)
    })

    output$p_umap <- renderPlot({
	p <- umap_plot()
	print(p)
    })

    # ---- PDF downloads ----
    output$download_pseudotime_pdf <- downloadHandler(
	filename = function() {
	    payload <- plot_payload()
	    req(payload)
	    sprintf("pseudotime_%s_%s.pdf",
		gsub("[^A-Za-z0-9]+", "_", payload$lineage),
		gsub("[^A-Za-z0-9]+", "_", payload$geneset_label))
	},
	content = function(file) {
	    # Use cairo_pdf for nicer text if available; otherwise use pdf()
	    # grDevices::cairo_pdf(file, width = 10, height = 6, onefile = TRUE)
	    grDevices::pdf(file, width = 10, height = 6, onefile = TRUE)
	    on.exit(grDevices::dev.off(), add = TRUE)
	    print(pseudotime_plot())
	}
    )

    output$download_umap_pdf <- downloadHandler(
	filename = function() {
	    payload <- plot_payload()
	    req(payload)
	    sprintf("umap_%s_%s.pdf",
		gsub("[^A-Za-z0-9]+", "_", payload$lineage),
		gsub("[^A-Za-z0-9]+", "_", payload$geneset_label))
	},
	content = function(file) {
	    # grDevices::cairo_pdf(file, width = 10, height = 8, onefile = TRUE)
	    grDevices::pdf(file, width = 10, height = 8, onefile = TRUE)
	    on.exit(grDevices::dev.off(), add = TRUE)
	    print(umap_plot())
	}
    )

    #---- raw table download ----
    output$download_raw_table <- downloadHandler(
	filename = function() {
	    payload <- plot_payload()
	    req(payload)

	    sprintf(
		"raw_table_%s_%s.tsv",
		gsub("[^A-Za-z0-9]+", "_", payload$lineage),
		gsub("[^A-Za-z0-9]+", "_", payload$geneset_label)
	    )
	},
	content = function(file) {
	    payload <- plot_payload()
	    req(payload)

	    d_out <- copy(payload$d)

	    # 可选：把 genes 和 missing genes 写成额外列，方便追踪
	    d_out[, genes_used := paste(payload$genes, collapse = ",")]
	    d_out[, missing_genes := paste(payload$missing, collapse = ",")]
	    d_out[, lineage := payload$lineage]
	    d_out[, geneset_label := payload$geneset_label]

	    fwrite(d_out, file, sep = "\t")
	}
    )



    output$p_pseudotime <- renderPlot({
	payload <- plot_payload()
	req(payload)

	d <- payload$d
	if (isTRUE(payload$facet)) {
	    ggplot(d, aes(x = pseudotime, y = gene_score)) +
		geom_point(alpha = 0.5, size = 1) +
		geom_smooth(method = "gam", se = FALSE, span = 0.05, n = 200) +
		facet_wrap(~marker, scales = "free_y") +
		theme_classic() +
		theme(axis.title = element_blank()) +
		ggtitle(sprintf("%s — %s pseudotime", payload$geneset_label, payload$lineage))
	} else {
	    plot_gene_vs_pseudotime(d, payload$geneset_label, payload$lineage, payload$selected_markers)
	}
    })

    output$p_umap <- renderPlot({
	payload <- plot_payload()
	req(payload)

	d <- payload$d[, .(UMAP1, UMAP2, gene_score, marker)] %>% na.omit()

	high_color = if (length(payload$selected_markers) == 1) {
	    marker_colors[payload$selected_markers]
	} else {
	    "red"
	}



	ggplot(d, aes(x = -UMAP1, y = UMAP2, color = gene_score)) +
	    geom_point(alpha = 0.5, size = 1) +
	    facet_wrap(~marker, scale = "free") +
	    scale_color_continuous(low = "lightgrey", high = high_color) +
	    theme_classic() +
	    theme(axis.title = element_blank()) +
	    ggtitle(sprintf("%s score on UMAP — %s lineage", payload$geneset_label, payload$lineage))
    })

    output$status <- renderText(status_msg())

    # nice error handling: show errors in status box instead of crashing plot
    observeEvent(input$submit, {
	tryCatch({
	    plot_payload()
	}, error = function(e) {
	    status_msg(paste("Error:", e$message))
	})
    })
}

shinyApp(ui, server, options = list(port = 4324, host = "23.121.124.95"))

