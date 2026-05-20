library(shiny)
library(plotly)
library(Seurat)
library(scales)

cellLasso <- function(rds, reduction = rds@misc$umap, dims = 1:2) {
  # Get rds data and save to data frame
  emb <- Embeddings(rds, reduction)[, dims]
  meta <- rds@meta.data
  ident <- Idents(rds)
  df <- data.frame(
    cell = rownames(emb),
    x = emb[,1],
    y = emb[,2],
    active_ident = ident,
    meta
  )
  
  # Plotly UI setup
  ui <- fluidPage(
    plotlyOutput("plot", height = "600px"),
    fluidRow(
      column(3, actionButton("save", "Save Selection")),
      column(3, actionButton("reset", "Reset Selection"))
    )
  )
  
  # Set umap colors to same as Seurat default
  default_cols <- hue_pal()(length(levels(ident)))
  
  server <- function(input, output, session) {
    
    output$plot <- renderPlotly({
      plot_ly(
        df,
        x = ~x,
        y = ~y,
        type = "scattergl",
        mode = "markers",
        color = ~ident,
        colors = default_cols,
        source = "lasso"
      ) %>%
        layout(dragmode = "lasso")
    })
    
    # Save button
    observeEvent(input$save, {
      sel <- event_data("plotly_selected", source = "lasso")
      
      if (is.null(sel)) {
        showNotification("No cells selected.", type = "error")
        return()
      }
      
      selected.cells <- df$cell[sel$pointNumber + 1]
      
      assign("selected.cells", selected.cells, envir = .GlobalEnv)
      
      showNotification(
        paste(length(selected.cells), "cells saved to selected.cells"),
        type = "message"
      )
    })
    
    # Reset button
    observeEvent(input$reset, {
      # Clear global variable selected.cells if it exists
      if (exists("selected.cells", envir = .GlobalEnv)) {
        rm(selected.cells, envir = .GlobalEnv)
      }
      
      # Reset the plot (re-render)
      output$plot <- renderPlotly({
        plot_ly(
          df,
          x = ~x,
          y = ~y,
          type = "scattergl",
          mode = "markers",
          source = "lasso"
        ) %>%
          layout(dragmode = "lasso")
      })
      
      showNotification("Selection reset.", type = "message")
    })
  }
  print("Selected cells saved to variable selected.cells")
  shinyApp(ui, server)
}

# MAIN ####
cellLasso(rds, reduction = "wnn.umap")


Idents(rds, cells= selected.cells) <- "Sample_01"
DimPlot(rds, reduction = "wnn.umap")
