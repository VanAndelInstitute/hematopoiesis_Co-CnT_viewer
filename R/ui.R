ui <- fluidPage(
    tags$head(
        tags$title("Trajectory Chromatin Explorer"),
        tags$style(HTML("
            :root {
                --bg: #f4f7fb;
                --panel: rgba(255,255,255,0.94);
                --panel-border: #d9e2ec;
                --ink: #102a43;
                --muted: #52606d;
                --accent: #0f6c8d;
                --accent-dark: #0b4f67;
                --warm: #d9a441;
                --shadow: 0 22px 60px rgba(15, 23, 42, 0.10);
            }
            body {
                background:
                    radial-gradient(circle at top right, rgba(15,108,141,0.18), transparent 28%),
                    linear-gradient(180deg, #eef4f8 0%, var(--bg) 38%, #f8fafc 100%);
                color: var(--ink);
                font-family: 'Trebuchet MS', Verdana, sans-serif;
            }
            .container-fluid {
                max-width: 1480px;
                padding: 28px 26px 44px;
            }
            .app-shell {
                display: flex;
                flex-direction: column;
                gap: 22px;
            }
            .hero-panel {
                padding: 28px 32px;
                border: 1px solid rgba(217,226,236,0.7);
                border-radius: 28px;
                background:
                    linear-gradient(135deg, rgba(255,255,255,0.96), rgba(240,244,248,0.92)),
                    linear-gradient(120deg, rgba(15,108,141,0.05), rgba(217,164,65,0.06));
                box-shadow: var(--shadow);
            }
            .hero-kicker {
                font-family: 'Trebuchet MS', Verdana, sans-serif;
                text-transform: uppercase;
                letter-spacing: 0.18em;
                font-size: 11px;
                font-weight: 700;
                color: var(--accent);
                margin-bottom: 10px;
            }
            .hero-title {
                font-size: 38px;
                line-height: 1.08;
                font-weight: 700;
                margin: 0 0 12px;
                color: var(--ink);
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .hero-copy {
                max-width: 920px;
                font-size: 16px;
                line-height: 1.7;
                color: var(--muted);
                margin: 0;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .content-grid {
                display: grid;
                grid-template-columns: 340px minmax(0, 1fr);
                gap: 22px;
                align-items: start;
            }
            .control-panel, .plot-card {
                background: var(--panel);
                border: 1px solid var(--panel-border);
                border-radius: 24px;
                box-shadow: var(--shadow);
            }
            .control-panel {
                padding: 24px 22px;
                position: sticky;
                top: 18px;
            }
            .panel-title {
                font-size: 22px;
                font-weight: 700;
                margin: 0 0 6px;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .panel-copy {
                font-size: 14px;
                line-height: 1.6;
                color: var(--muted);
                margin: 0 0 18px;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .form-group > label {
                font-family: 'Trebuchet MS', Verdana, sans-serif;
                font-size: 12px;
                text-transform: uppercase;
                letter-spacing: 0.12em;
                color: var(--muted);
                margin-bottom: 8px;
            }
            .form-control, .selectize-input, .selectize-dropdown, textarea {
                border-radius: 14px !important;
                border-color: #cbd5e1 !important;
                box-shadow: none !important;
                font-size: 14px;
            }
            .form-control:focus, .selectize-input.focus {
                border-color: rgba(15,108,141,0.65) !important;
                box-shadow: 0 0 0 3px rgba(15,108,141,0.12) !important;
            }
            .checkbox label {
                color: var(--ink);
                font-size: 14px;
            }
            .btn-primary {
                width: 100%;
                border: 0;
                border-radius: 999px;
                background: linear-gradient(135deg, var(--accent) 0%, var(--accent-dark) 100%);
                padding: 13px 18px;
                font-family: 'Trebuchet MS', Verdana, sans-serif;
                font-size: 13px;
                font-weight: 700;
                text-transform: uppercase;
                letter-spacing: 0.12em;
            }
            .btn-default {
                border-radius: 999px;
                border: 1px solid #bcccdc;
                background: white;
                color: var(--ink);
                font-family: 'Trebuchet MS', Verdana, sans-serif;
                font-size: 12px;
                font-weight: 700;
                text-transform: uppercase;
                letter-spacing: 0.08em;
                padding: 10px 16px;
            }
            .status-box {
                margin-top: 18px;
                padding: 16px 18px;
                border-radius: 16px;
                background: linear-gradient(180deg, #f8fbfd 0%, #f1f5f9 100%);
                border: 1px solid #d9e2ec;
                font-size: 13px;
                line-height: 1.6;
                color: var(--ink);
                white-space: pre-wrap;
                font-family: 'Trebuchet MS', Verdana, sans-serif;
            }
            .results-grid {
                display: grid;
                grid-template-columns: 1fr;
                gap: 22px;
            }
            .plot-card {
                padding: 22px 24px 18px;
            }
            .dataset-details {
                background: var(--panel);
                border: 1px solid var(--panel-border);
                border-radius: 24px;
                box-shadow: var(--shadow);
                padding: 0;
                overflow: hidden;
            }
            .dataset-details summary {
                list-style: none;
                cursor: pointer;
                padding: 20px 24px;
                display: flex;
                align-items: center;
                justify-content: space-between;
                gap: 16px;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
                font-size: 22px;
                font-weight: 700;
                color: var(--ink);
            }
            .dataset-details summary::-webkit-details-marker {
                display: none;
            }
            .dataset-details summary::after {
                content: '+';
                font-family: 'Trebuchet MS', Verdana, sans-serif;
                font-size: 26px;
                line-height: 1;
                color: var(--accent-dark);
            }
            .dataset-details[open] summary::after {
                content: '-';
            }
            .dataset-copy {
                padding: 0 24px 22px;
                color: var(--muted);
                font-size: 14px;
                line-height: 1.7;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .dataset-copy p:last-child {
                margin-bottom: 0;
            }
            .dataset-copy a {
                color: var(--accent-dark);
            }
            .dataset-figure {
                margin: 18px 0 20px;
            }
            .dataset-figure img {
                display: block;
                width: 100%;
                max-width: 500px;
                height: auto;
                margin: 0 auto;
                border-radius: 18px;
                border: 1px solid var(--panel-border);
                background: #fff;
            }
            .plot-card-header {
                display: flex;
                justify-content: space-between;
                align-items: flex-end;
                gap: 16px;
                margin-bottom: 14px;
            }
            .plot-title {
                font-size: 23px;
                font-weight: 700;
                margin: 0;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .plot-copy {
                margin: 4px 0 0;
                color: var(--muted);
                font-size: 14px;
                line-height: 1.6;
                max-width: 760px;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
            }
            .download-row {
                display: flex;
                flex-wrap: wrap;
                gap: 10px;
            }
            .shiny-plot-output {
                border-radius: 18px;
                overflow: hidden;
            }
            .app-footer {
                padding: 4px 6px 0;
                text-align: center;
                color: var(--muted);
                font-size: 13px;
                line-height: 1.7;
                font-family: 'Trebuchet MS', Verdana, sans-serif;
            }
            .app-footer a {
                color: var(--accent-dark);
            }
            @media (max-width: 1100px) {
                .content-grid {
                    grid-template-columns: 1fr;
                }
                .control-panel {
                    position: static;
                }
                .plot-card-header {
                    flex-direction: column;
                    align-items: flex-start;
                }
            }
        "))
    ),
    div(
        class = "app-shell",
        div(
            class = "hero-panel",
            div(class = "hero-kicker", "Hematopoiesis Histone modifications dynamics Viewer"),
            h1(class = "hero-title", "Chromatin Dynamics Across Hematopoietic Lineages"),
            p(
                class = "hero-copy",
                "Single cell CoCut&Tag Bone marrow Hematopoietic stem cell histone modification browser. "
            )
        ),
        div(
            class = "content-grid",
            div(
                class = "control-panel",
                h2(class = "panel-title", "Analysis Design"),
                p(
                    class = "panel-copy",
                    "Define a lineage, a gene or a gene panel. The app computes averaged gene scores and renders synchronized pseudotime and UMAP embedding views."
                ),
                selectInput(
                    "lineage", "Lineage / trajectory",
                    choices = names(trajectory_map),
                    selected = "B cell"
                ),
                textAreaInput(
                    "genes", "Genes (one per line)",
                    value = "PAX5",
                    placeholder = "PAX5",
                    rows = 7,
                    width = "100%"
                ),
                helpText(sprintf("Maximum %d genes per query.", MAX_GENES_PER_QUERY)),
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
                            "min-max (0-1)" = "minmax"
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
                actionButton("submit", "Update Figures", class = "btn-primary"),
                div(class = "status-box", textOutput("status"))
            ),
            div(
                class = "results-grid",
                div(
                    class = "plot-card",
                    div(
                        class = "plot-card-header",
                        div(
                            h3(class = "plot-title", "Pseudotime Profile"),
                            p(class = "plot-copy", "Histone modification gene score versus lineage pseudotime.")
                        ),
                        div(
                            class = "download-row",
                            downloadButton("download_pseudotime_pdf", "Pseudotime PDF"),
                            downloadButton("download_raw_table", "Raw Table TSV")
                        )
                    ),
                    plotOutput("p_pseudotime", height = "760px")
                ),
                div(
                    class = "plot-card",
                    div(
                        class = "plot-card-header",
                        div(
                            h3(class = "plot-title", "Embedding View"),
                            p(class = "plot-copy", "UMAP projection colored by histone modification score score.")
                        ),
                        div(
                            class = "download-row",
                            downloadButton("download_umap_pdf", "UMAP PDF")
                        )
                    ),
                    plotOutput("p_umap", height = "760px")
                ),
                tags$details(
                    class = "dataset-details",
                    tags$summary("About the dataset"),
                    div(
                        class = "dataset-copy",
                        p("Single-cell CoCut&Tag bone marrow hematopoietic stem cell histone modification browser."),
                        div(
                            class = "dataset-figure",
                            tags$img(
                                src = "data/fig_dataset_diagram.png",
                                alt = "Dataset diagram for the hematopoietic histone modification browser"
                            )
                        ),
                        p(
                            "TBD: the dataset details, including experimental design, data processing, and analysis methods will be described here. For now, please refer to the associated publication (cite the preprint) or contact the lab for more information."
                        )
                    )
                )
            )
        ),
        div(
            class = "app-footer",
            "Janssens Lab, ",
            tags$a(href = "mailto:Derek.Janssens@vai.org", "Derek.Janssens@vai.org"),
            " | ",
            tags$a(
                href = "https://janssenslab.vai.org/",
                target = "_blank",
                rel = "noopener noreferrer",
                "https://janssenslab.vai.org/"
            )
        )
    )
)
