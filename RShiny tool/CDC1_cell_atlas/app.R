#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(shinythemes)
library("Seurat")
library("plotly")

SeuratObjT <- readRDS("C:/Users/samue/Desktop/Stage VIB/Pre and Post processing Results/SAM2/Post-process/SeuratObj_Post-Process_CDC1_SAM2.rds")
SeuratObjT <- DietSeurat(SeuratObjT,
              counts = FALSE,
              data = TRUE,
              scale.data = FALSE,
              assays = "RNA",
              dimreducs = c("RNA_umap")
)
# List of Genes
#genes <- rownames(seuratObjT)
genes <- c("Test","Ccr7")

# Define UI for application that draws a histogram

############## Base of UI #############################
ui <-  page_navbar(theme = shinytheme("superhero"),
  title = img(src="irc_logo_transparant.png",  
              height = "30px"),
  bg = "#8EE5EE",
  inverse = TRUE,

############ CCS style #######################
tags$head(
  tags$style(HTML("
            .Title{
              font-size: 32px;
              margin-bottom: 5px;
            }
            .body {
              font-family: sans-serif;
            }
         .card{
  position: relative;
  width: 100%;
  max-width: 800px;
  height: 700px;

  margin: 2em auto;
  padding: 1em;

  overflow: hidden;

  display: flex;
  flex-direction: column;
}

.card .plotly {
  flex: 1;
  min-height: 0;
}
            .input{
              position:relative;
              padding: 2em;
              margin: 2em auto;
            
            }
            
                    "))
),

################## Home Page ##################################
  nav_panel(title = "Home",
      tags$div(
        class = "Title",
        "Home"
      ),
      
      tags$div( class= "card",
          p("This is cell atlas for CDC1 based on these Datasets:"),
          
      p("CITEseq_Test: SAM2 and SAM3",tags$br(),
      "CITEseq_Final: SAM05 and SAM06",tags$br(),
      "CITEseq_Notch: SAM016\n",       tags$br(),
      "CITEseq_LNP: VBO004-VBO012\n",  tags$br(),
      "CITEseq_Toxo: JVE008 JVE010.") )
  ),
 
##################################################################
############## Feature Plot Page ################################
  nav_panel(title = "Feature plots",
            
            p("Select gene for feature plot."),
            
            tags$div( class = "input",
                      
            selectizeInput(
              inputId = "gene",
              label = "Select or type a gene",
              choices = genes,
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE   # allows typing custom values
              )),
            verbatimTextOutput("selected_gene")
            ),
            
            tags$div(
              plotlyOutput("feature_plot", height = "700px")
            )
          ),
  
#############################################################################
################# TBD ######################################################
  nav_panel(title = "Three", p("Third page content.")),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Posit", href = "https://posit.co")),
    nav_item(tags$a("Shiny", href = "https://shiny.posit.co"))
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {

  output$selected_gene <- renderPrint({input$gene})
  
  output$feature_plot <- renderPlotly({ 
    
    req(input$gene)
    
    umap <- SeuratObjT@reductions$RNA_umap@cell.embeddings
    
    expr <- FetchData(SeuratObjT, vars = input$gene)
    
    df <- data.frame(
      UMAP_1 = umap[, "RNAumap_1"],
      UMAP_2 = umap[, "RNAumap_2"],
      expr = expr[, 1]
    )
    
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr,
      type = "scatter",
      mode = "markers",
      marker = list(size = 3)
    )
  })
}
# Run the application 
shinyApp(ui = ui, server = server)
