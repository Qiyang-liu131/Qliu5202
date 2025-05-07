library(shiny)  # Loads the Shiny library, which is used to build interactive web applications in R.
library(Seurat) # Loads the Seurat library for single-cell RNA-seq data analysis.
library(ggplot2) # Loads ggplot2, a library for creating visualizations in R.

# Raw GitHub URL to the .rds file
url <- "https://raw.githubusercontent.com/Qiyang-liu131/Qliu5202/main/merged_KOmouse_seurat.rds"
# Load the RDS file directly

master_seurat <- readRDS(url(url)) 

# User Interface (UI) definition for the Shiny application.
ui <- fluidPage(
  titlePanel("Spermatogenesis Data UMAP"), # Adds a title to the application.
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput(
        "selectedStages", # Input ID.
        "Select Stages to Display:", # Label for the input.
        choices = unique(master_seurat$stage), # Provides a list of unique stages from the Seurat object for user selection.
        selected = unique(master_seurat$stage) # By default, selects all stages.
      ),
      actionButton("updatePlot", "Update Plot") # Adds a button to update the UMAP plot based on user input.
    ),
    mainPanel(
      plotOutput("umapPlot") # Placeholder for the UMAP plot output.
    )
  )
)

# Server logic for the Shiny application.
server <- function(input, output) {
  
  # Reactive expression to filter the Seurat object based on user-selected stages.
  filtered_seurat <- reactive({
    req(input$selectedStages) # Ensures the input is available before proceeding.
    subset(master_seurat, subset = stage %in% input$selectedStages) # Filters the Seurat object to include only selected stages.
  })
  
  # Renders the UMAP plot based on the filtered Seurat object.
  output$umapPlot <- renderPlot({
    input$updatePlot # Triggers the plot update when the action button is clicked.
    isolate({ # Ensures that only the action button triggers the plot update, not the reactive input changes.
      DimPlot(filtered_seurat(), reduction = "umap", group.by = "stage") + # Generates a UMAP plot grouped by stage.
        ggtitle("UMAP by Selected Stages") # Adds a title to the plot.
    })
  })
}

# Combines the UI and server to create and launch the Shiny application.
shinyApp(ui = ui, server = server)
