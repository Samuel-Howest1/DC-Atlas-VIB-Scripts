#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library("Seurat")
library("plotly")
library(bslib)
library(bootstrap)
library(shinyWidgets)
library(DT)
library(shinycssloaders)
################### Loading Objects #######################################
SeuratObjT <- readRDS("C:/Users/irc/Desktop/Harmony_Treat_Exp_VIS.rds")
SeuratObjT <- JoinLayers(SeuratObjT,assay = "RNA")
Idents(SeuratObjT) <- SeuratObjT$sctype_classification
################################################################################

################### Creating List of of selection for plots #########
# genes <- rownames(SeuratObjT)
genes <- c("Test","Ccr7","Cd81")
condition <- c(
  paste("Orig.ident:",unique(SeuratObjT@meta.data$orig.ident)),
  paste("Treatment:",unique(SeuratObjT@meta.data$treatment)),
  paste("Experiment:",unique(SeuratObjT@meta.data$experiment))
               )
metadata <- c("orig.ident","seurat_clusters","sctype_classification", "scDblFinder_class")
################################################################################
##################### Spliting up the data for improved speed #################
vis_data <- list(
  umap = SeuratObjT@reductions$harmony_umap@cell.embeddings,
  meta = SeuratObjT@meta.data,
  expr = GetAssayData(SeuratObjT,assay = "RNA",layer = "data"))

umap <- vis_data$umap
meta <- vis_data$meta
expr <- vis_data$expr

##############################################################################
########################### Info Table | Home Page ###########################
# X-axis (columns)
report_cols <-c(
    "JVE008", "JVE010", "SAM016",
    "SAM05", "SAM06", "SAM2",
    "SAM3", "VBO004", "VBO005",
    "VBO006", "VBO007", "VBO008",
    "VBO009", "VBO010", "VBO011",
    "VBO012"
  )
experiment_map <- c(
  "CITEseq_Test",
  "CITEseq_Test",
  "CITEseq_Final",
  "CITEseq_Final",
  "CITEseq_Notch",
  "CITEseq_LNP_WT",
  "CITEseq_LNP_eLNPs",
  "CITEseq_LNP_pIC_LNPs",
  "CITEseq_LNP_CpG",
  "CITEseq_LNP_pIC",
  "CITEseq_LNP_eLNPs",
  "CITEseq_LNP_pIC_LNPs",
  "CITEseq_LNP_CpG",
  "CITEseq_LNP_pIC",
  "CITEseq_Toxo",
  "CITEseq_Toxo"
)
# Y-axis (rows)
conditions <- c(as.character(
  tags$div(
    tags$img(src = "Orig_idents.png", height = "30px")," Orig.idents"
    )),
  as.character(
    tags$div(tags$img(src = "mouse_icon.png", height = "30px"),
      " Treatment")),
  as.character(
    tags$div(
      tags$img(src = "lab.png", height = "30px")," Experiment"
    ))
)
tbl <- data.frame(
  Condition = conditions,
  matrix(
    "",
    nrow = length(conditions),
    ncol = length(report_cols),
    dimnames = list(NULL, report_cols)
  ),
  check.names = FALSE)
##############################################################################
############## Base of UI ###################################################
ui <-  page_navbar(theme = shinytheme("united"),
  title = img(src="irc_logo_transparant.png",  
              height = "30px"),
  bg = "#8EE5EE",
  inverse = TRUE,

############ CCS style ###################################################
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
              overflow: vissible;
              display: flex;
              flex-direction: column;
         }
         
        .table{
              width: 1200px;
              margin-left:0px;
        }

        .card .plotly {
              flex: 1;
              min-height: 0;
        }

        .input{
            position:relative;
            width: 100%;
            padding: 2em;
            margin: 2em auto;
            
            }
            
                    "))
),

################## Home Page ##################################
  nav_panel(title = "Home",
      tags$div(
        class = "Title",
        h1("Home")
      ),
      
      tags$div(
          p("This is cell atlas for CDC1 based on these Datasets:"),
          
      p("CITEseq_Test: SAM2 and SAM3",tags$br(),
      "CITEseq_Final: SAM05 and SAM06",tags$br(),
      "CITEseq_Notch: SAM016\n",       tags$br(),
      "CITEseq_LNP: VBO004-VBO012\n",  tags$br(),
      "CITEseq_Toxo: JVE008 JVE010.")),
      
      div(class ="table",
        tableOutput("report_table")
      )
  ),
 
##################################################################
############## Feature Plot Page ################################
  nav_panel(title = "Gene plots",
# ------------------ Selectors -----------------------------------
card(
          max_height = "150px",
            p("Select gene for feature plot."),
layout_columns(              
            selectizeInput(
              inputId = "gene",
              label = "Select or type a gene",
              choices = genes,
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE,   # allows typing custom values
                dropdownParent = "body"
              )),
            selectizeInput(
              inputId = "cond",
              label = "Select or type a gene",
              choices = condition,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE,   # allows typing custom values
                dropdownParent = "body"
              ),)

            # radioButtons(
            #   inputId = "cond",
            #   label = "Condition",
            #   choices = c(
            #     "Orig.idebts" = "orig",
            #     "Treatment" = "treat",
            #     "WT Only" = "wt"
            #   ),
            #   selected = "orig",
            #   inline = TRUE
            # )
          ,col_widths = c(6, 6))
),

# ------------------------------------------------------------------

  layout_columns(
          
        card(  
          withSpinner(
          plotlyOutput("feature_plot", height = "700px",width = "700px")
          )),
        card(
            card_header("Violin Plot"),
            plotOutput("ViolinPlot", height = "700px"),
            
            plotOutput("DotPlot", height = "700px")
            
          )
          ),col_widths = c(8, 4)
        ),
  
#############################################################################
################# Cell Metdata ######################################################
  nav_panel(title = "Cell Metadata", 
            p("Select gene for feature plot."),
            tags$div( class = "input",
                      
            selectizeInput(
              inputId = "meta",
              label = "Select Label for metadata",
              choices = genes,
              selected = NULL,
              multiple = F,
              options = list(
              placeholder = 'Type a gene...',
              create = TRUE   # allows typing custom values
              )),
            
            selectizeInput(
              inputId = "cond",
              label = "Select or type a gene",
              choices = condition,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE,   # allows typing custom values
                dropdownParent = "body"
              ),)

            ),
            
            withSpinner(
            plotlyOutput("metaplot", height = "1000px",width = "1000px")
            )
            # 
            # tags$div(
            #   plotOutput("Dimplot", height = "700px",width = "700px")
            # )
  ),


###############################################################################
############################ Gene Comparison #################################

nav_panel(title = "Gene Comparison", 
          p("Select gene for feature plot."),
          tags$div( class = "input",
                    
                    selectizeInput(
                      inputId = "gene_1",
                      label = "Select Label for genes",
                      choices = genes,
                      selected = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = 'Type a gene...',
                        create = TRUE   # allows typing custom values
                      )),
                    
                    selectizeInput(
                      inputId = "gene_2",
                      label = "Select Label for genes",
                      choices = genes,
                      selected = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = 'Type a gene...',
                        create = TRUE   # allows typing custom values
                      )),
          ),
           
          tags$div(
            plotOutput("Comparison", height = "700px",width = "700px")
          )
),

##############################################################################
############################## Contact #######################################
nav_panel(title = "Contact", 
          p("Select gene for feature plot."),
          tags$div( class = "input",
          ),
          
          tags$div(
          )
),

################################################################################
################################ Links ########################################
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Posit", href = "https://posit.co")),
    nav_item(tags$a("Shiny", href = "https://shiny.posit.co"))
  )
)

#####################################################################################
############################      SERVER     ########################################
#####################################################################################
# Define server logic required to draw a histogram
server <- function(input, output) {
  

######################################################################################
  #Home Table
  output$report_table <- renderTable(
    tbl,
    striped = TRUE,
    bordered = TRUE,
    spacing = "s",
    # Makes it so that text within the table will also be treated as a html elements
    sanitize.text.function = function(x) x
  )
  
######################################################################################
  output$selected_gene <- renderPrint({input$gene})
  
  output$feature_plot <- renderPlotly({ 
    
    req(input$gene)
    
    # Cell coordinetes (already created)
    # umap <- SeuratObjT@reductions$harmony_umap@cell.embeddings
    
    # Gene 
    expr_val <-  as.numeric(expr[input$gene, ])
    # Celltype
    cell <- meta$sctype_classification
    
    # Table for plotting
    df <- data.frame(
      UMAP_1 = umap[, "harmonyumap_1"],
      UMAP_2 = umap[, "harmonyumap_2"],
      expr = expr_val
    )
    df$Celltype <- meta$sctype_classification
    # To male Hover Text of cells
    df$hover <- paste(
      rownames(df),
      "<br>Expression:",
      round(df$expr, 2),
      "<br>Celltype:",df$Celltype
    )
    # Making a Plotly plot that resembles a featueplot
    plot_ly(
      df,
      x = ~UMAP_1,
      y = ~UMAP_2,
      color = ~expr,
      type = "scattergl",
      mode = "markers",
      text = ~hover,
      #Removing Coordintates
      hoverinfo= "text",
      marker = list(size = 3)
    )
  })
  
  output$ViolinPlot <- renderPlot({
    
    req(input$gene)
    
    VlnPlot(
      SeuratObjT,
      features = input$gene,
      pt.size = 0
    )
  })
  
  # Info verbeter zodat het duielijkt is 
  output$DotPlot <- renderPlot({
    DotPlot(
      SeuratObjT,
      features = genes,
      group.by = "sctype_classification"
    )
  })
  
#################################################################################

  ########## Cell Metadata Start button ############   
  PlotNameList <-eventReactive(input$start,{
    Count <- 1
    for (i in input$cond) {
    word <- strsplit(i,":")[[1]]
    
    }
  })
  
-------------------------------------------------------------------------
 ########## Creating the Cell Metadata Dynamic plots ################
  output$metaplot <- renderPlotly({
      
      req(PlotNameList())
      req(input$meta)
      expr_val <-  as.numeric(expr[input$meta, ])
      # Celltype
      cell <- meta$sctype_classification
      # Table for plotting
      df <- data.frame(
        UMAP_1 = umap[, "harmonyumap_1"],
        UMAP_2 = umap[, "harmonyumap_2"],
        expr = expr_val
      )
      
      df$Exp <- meta$experiment
      df$WT <- meta$treatment
      Count <- 1
      
      for (i in PlotNameList){
        df_Filter <- df[df == i,]
        p <-plot_ly(
            df_Filter,
            x = ~UMAP_1,
            y = ~UMAP_2,
            color = ~expr,
            type = "scattergl",
            mode = "markers",
            marker = list(size = 3)
          )
        PlotsList[[Count]] <- p
        Count <- Count + 1

      }
      
      # df_Filter_WT <- df[df$WT == "WT",]
      # p1 <-plot_ly(
      #   df_Filter_WT,
      #   x = ~UMAP_1,
      #   y = ~UMAP_2,
      #   color = ~expr,
      #   type = "scattergl",
      #   mode = "markers",
      #   marker = list(size = 3)
      # )
      # 
      # df_Filter_Test <- df[df$WT == "Test",]
      # p2 <- plot_ly(
      #   df_Filter_Test,
      #   x = ~UMAP_1,
      #   y = ~UMAP_2,
      #   color = ~expr,
      #   type = "scattergl",
      #   mode = "markers",
      #   marker = list(size = 3)
      # )
      # 
      # df_Filter_JVE <- df[df$Exp == "CITEseq_Toxo",]
      # p3 <- plot_ly(
      #   df_Filter_JVE,
      #   x = ~UMAP_1,
      #   y = ~UMAP_2,
      #   color = ~expr,
      #   type = "scattergl",
      #   mode = "markers",
      #   marker = list(size = 3)
      # )
      # subplot(
      #   p1,
      #   p2,
      #   p3,
      #   nrows = 1
      # )
    })
  
  
#################################################################################
# ADD COUNT OF HOW MANY CELL PER GENE
  output$Comparison <- renderPlot({
    FeatureScatter(
      SeuratObjT,
      feature1 = input$gene_1,
      feature2 = input$gene_2
    )
  })
  
}
# Run the application 
shinyApp(ui = ui, server = server)


