#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(shiny)
    library(data.table)
    library(ggplot2)
    library(magrittr)
})

for (script in c("helpers.R", "globals.R", "ui.R", "server.R")) {
    source(file.path("R", script), local = FALSE)
}

addResourcePath("data", "data")

shinyApp(ui, server, options = list(port = 4324, host = "23.121.124.95"))
