server <- function(input, output, session) {
    status_msg <- reactiveVal("Choose gene/lineage/markers, then click Submit / Plot.")
    query_cache <- make_lru_cache(QUERY_CACHE_MAX)

    plot_payload <- eventReactive(input$submit, {
        lineage <- isolate(input$lineage)
        normalize <- isolate(input$normalize_genes)
        norm_method <- isolate(input$norm_method)
        genes_txt <- isolate(input$genes)
        allow_missing_genes <- isolate(input$allow_missing_genes)
        markers <- isolate(input$markers)
        facet <- isolate(input$facet_markers)

        genes <- parse_genes_from_text(genes_txt)
        selected_markers <- intersect(markers, available_markers)

        if (length(genes) == 0) {
            stop("Please input at least one gene (one per line).")
        }
        if (length(genes) > MAX_GENES_PER_QUERY) {
            stop(sprintf(
                "Please limit the query to %d genes or fewer. You submitted %d genes.",
                MAX_GENES_PER_QUERY,
                length(genes)
            ))
        }
        if (length(selected_markers) == 0) {
            stop("Please select at least one marker.")
        }

        cache_key <- paste(
            lineage,
            paste(sort(genes), collapse = ","),
            paste(sort(selected_markers), collapse = ","),
            normalize,
            norm_method,
            allow_missing_genes,
            sep = "|"
        )

        if (query_cache$exists(cache_key)) {
            payload <- query_cache$get(cache_key)
            payload$facet <- facet
            status_msg(make_status_message(
                "OK (cached)",
                payload$geneset_label,
                lineage,
                length(unique(payload$d$marker)),
                nrow(payload$d)
            ))
            return(payload)
        }

        pt <- make_pseudotime_dt(lineage)
        gs_obj <- make_geneset_score_dt(
            genes,
            selected_markers,
            allow_missing = allow_missing_genes,
            normalize = normalize,
            norm_method = norm_method
        )
        gs <- gs_obj$dt

        d <- merge(gs, pt, by = "cell_id")
        d <- merge(d, umap_k27, by = "cell_id", all.y = TRUE)
        if (nrow(d) == 0) {
            stop("No data after merging gene scores with pseudotime.")
        }

        geneset_label <- sprintf("Avg(%d genes)", length(gs_obj$present))
        if (length(gs_obj$missing) > 0) {
            geneset_label <- sprintf("%s | missing:%d", geneset_label, length(gs_obj$missing))
        }

        status_text <- make_status_message(
            "OK",
            geneset_label,
            lineage,
            length(unique(d$marker)),
            nrow(d)
        )
        if (length(gs_obj$missing) > 0) {
            status_text <- paste0(
                status_text,
                "\nUsed: ",
                paste(gs_obj$present, collapse = ", "),
                "\nMissing: ",
                paste(gs_obj$missing, collapse = ", ")
            )
        }
        status_msg(status_text)

        payload <- list(
            d = d,
            genes = gs_obj$present,
            missing = gs_obj$missing,
            geneset_label = geneset_label,
            lineage = lineage,
            selected_markers = selected_markers,
            facet = facet
        )

        query_cache$set(cache_key, payload)
        payload
    }, ignoreInit = TRUE)

    pseudotime_plot <- reactive({
        payload <- plot_payload()
        req(payload)
        d <- payload$d

        if (isTRUE(payload$facet)) {
            ggplot(d, aes(x = pseudotime, y = gene_score)) +
                geom_point(alpha = 0.35, size = 0.85, color = "#486581") +
                geom_smooth(method = "gam", se = FALSE, color = "#0f6c8d", linewidth = 1.05, span = 0.08, n = 300) +
                facet_wrap(~marker, scales = "free_y", ncol = 2) +
                labs(
                    title = "Gene-Set Dynamics Across Pseudotime",
                    subtitle = plot_subtitle(payload$lineage, payload$geneset_label, length(payload$selected_markers)),
                    x = "Pseudotime",
                    y = "Gene Set Score"
                ) +
                scientific_theme()
        } else {
            plot_gene_vs_pseudotime(
                d,
                payload$geneset_label,
                payload$lineage,
                payload$selected_markers
            )
        }
    })

    umap_plot <- reactive({
        payload <- plot_payload()
        req(payload)

        d <- payload$d[, .(UMAP1, UMAP2, gene_score, marker)] %>% na.omit()

        high_color <- if (length(payload$selected_markers) == 1) {
            marker_colors[payload$selected_markers]
        } else {
            "#0f6c8d"
        }

        ggplot(d, aes(x = -UMAP1, y = UMAP2, color = gene_score)) +
            geom_point(alpha = 0.72, size = 0.95) +
            facet_wrap(~marker, scales = "fixed", ncol = 2) +
            scale_color_gradient(low = "#d9e2ec", high = high_color) +
            labs(
                title = "UMAP Distribution of Gene-Set Scores",
                subtitle = plot_subtitle(payload$lineage, payload$geneset_label, length(payload$selected_markers)),
                x = "UMAP 1",
                y = "UMAP 2",
                color = "Score"
            ) +
            coord_equal() +
            scientific_theme()
    })

    output$p_pseudotime <- renderPlot({
        print(pseudotime_plot())
    }, res = 120)

    output$p_umap <- renderPlot({
        print(umap_plot())
    }, res = 120)

    output$download_pseudotime_pdf <- downloadHandler(
        filename = function() {
            payload <- plot_payload()
            req(payload)
            make_export_filename("pseudotime", payload, "pdf")
        },
        content = function(file) {
            grDevices::pdf(file, width = 12.5, height = 8.2, onefile = TRUE)
            on.exit(grDevices::dev.off(), add = TRUE)
            print(pseudotime_plot())
        }
    )

    output$download_umap_pdf <- downloadHandler(
        filename = function() {
            payload <- plot_payload()
            req(payload)
            make_export_filename("umap", payload, "pdf")
        },
        content = function(file) {
            grDevices::pdf(file, width = 12.5, height = 8.4, onefile = TRUE)
            on.exit(grDevices::dev.off(), add = TRUE)
            print(umap_plot())
        }
    )

    output$download_raw_table <- downloadHandler(
        filename = function() {
            payload <- plot_payload()
            req(payload)
            make_export_filename("raw_table", payload, "tsv")
        },
        content = function(file) {
            payload <- plot_payload()
            req(payload)

            d_out <- copy(payload$d)
            d_out[, genes_used := paste(payload$genes, collapse = ",")]
            d_out[, missing_genes := paste(payload$missing, collapse = ",")]
            d_out[, lineage := payload$lineage]
            d_out[, geneset_label := payload$geneset_label]

            fwrite(d_out, file, sep = "\t")
        }
    )

    output$status <- renderText(status_msg())

    observeEvent(input$submit, {
        tryCatch({
            plot_payload()
        }, error = function(e) {
            status_msg(paste("Error:", e$message))
        })
    })
}
