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
library(ggplot2)
library(shinycssloaders)
################### Loading Objects #######################################
SeuratObjT <- readRDS("C:/Users/irc/Desktop/Harmony_Treat_Exp_VIS.rds")
SeuratObjT <- JoinLayers(SeuratObjT,assay = "RNA")
Idents(SeuratObjT) <- SeuratObjT$sctype_classification
################################################################################
##################### Spliting up the data for improved speed #################
vis_data <- list(
  umap = SeuratObjT@reductions$harmony_umap@cell.embeddings,
  meta = SeuratObjT@meta.data,
  expr = GetAssayData(SeuratObjT,assay = "RNA",layer = "data"))

umap <- vis_data$umap
meta <- vis_data$meta
expr <- vis_data$expr

### Removing Seuratobject
rm(SeuratObjT)
rm(vis_data)
gc()
################### Creating List of of selection for plots #########
condition <- c(
  paste("Orig.ident:",unique(meta$orig.ident)),
  paste("Treatment:",unique(meta$treatment)),
  paste("Experiment:",unique(meta$experiment))
)
metadata <- c("orig.ident","seurat_clusters","sctype_classification", "scDblFinder_class")
################################################################################

# List of Genes
genes <- rownames(expr)
#genes <- c("Test","Ccr7","Cd81")

########### Avg and Percentage expression of Genes for Dotplot later ##########
# avg <- AverageExpression(SeuratObjT,group.by = "sctype_classification",features = genes,layer = "data")$RNA
# 
# pct_group <- meta$sctype_classification
# pct_exp <- function(gene, expr, pct_group) {
#   tapply(expr[gene, ] > 0, pct_group, mean) * 100
#   }
# pct_list <- lapply(genes, function(g) pct_exp(g, expr, group))
# 
# pct_mat <- do.call(cbind, pct_list)
# colnames(pct_mat) <- genes
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
          .metaplots{
          margin:50px;
          padding; 50px;
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
layout_columns(              
            selectizeInput(
              inputId = "gene",
              label = "Select gene for feature plot",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE,   # allows typing custom values
                dropdownParent = "body"
              ))
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
            
            selectizeInput(
              inputId = "Dimplot",
              choices = metadata,
              label = "Select MetaData for Dimplot",
              selected = NULL,
              multiple = FALSE,
              options = list(
                placeholder = 'Type a gene...',
                create = TRUE,   # allows typing custom values
                dropdownParent = "body"
              )),
            
            plotOutput("DimPlotMeta", height = "700px")
            
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
              choices = NULL,
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
            radioButtons(
              inputId = "colour_by",
              label = "Colour by",
              choices = c(
                "Gene expression" = "gene",
                "Cell type" = "celltype",
                "Cluster" = "cluster",
                "Treatment" = "treatment",
                "Experiment" = "experiment",
                "Orig.ident" = "orig.ident"
              ),
              selected = "gene",
              inline = TRUE
            ),
            actionButton("start", "Generate plots"),
            
            withSpinner(
            plotlyOutput("metaplot", height = "700px",width = "1200px")
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
                      choices = NULL,
                      selected = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = 'Type a gene...',
                        create = TRUE   # allows typing custom values
                      )),
                    
                    selectizeInput(
                      inputId = "gene_2",
                      label = "Select Label for genes",
                      choices = NULL,
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
server <- function(input, output,session) {
  
  
#################### LIST OF GENE ###################################################
  updateSelectizeInput(session,"gene",choices = genes,server = TRUE)
  
  updateSelectizeInput(session,"meta",choices = genes,server = TRUE)
  
  updateSelectizeInput(session,"gene_1",choices = genes,server = TRUE)
  
  updateSelectizeInput(session,"gene_2",choices = genes,server = TRUE)

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
    
    # Umap 
    df <- data.frame(
      UMAP_1 = umap[, "harmonyumap_1"],
      UMAP_2 = umap[, "harmonyumap_2"]
    )
    
    df$expr <- as.numeric(expr[input$gene, ])
    
    # Table for plotting
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
    
    # VlnPlot(
    #   SeuratObjT,
    #   features = input$gene,
    #   pt.size = 0
    # )
    
    df_violin <- data.frame(
      expr= as.numeric(expr[input$gene,]),
      celltype = meta$sctype_classification
    )
    
    ggplot(df_violin, aes(x = celltype, y= expr)) +
      geom_violin()+
      theme_bw()
  })
  
  # Info verbeter zodat het duielijkt is 
  # output$DotPlot <- renderPlot({
  #   DotPlot(
  #     SeuratObjT,
  #     features = genes,
  #     group.by = "sctype_classification"
  #   )
  
  output$DimplotMeta <- renderPlot({
    
  ggplot(df,aes(x = x, y = y, color= input$Dimplot))+
      geom_point(size= 0.5)+
      theme_classic()+
      labs(x= "UMAP_1", y= "UMAP_2")
  })
  
#################################################################################

  ########## Cell Metadata Start button ############   
  PlotNameList <-eventReactive(input$start,{
      input$cond
  })
  
# -------------------------------------------------------------------------
 ########## Creating the Cell Metadata Dynamic plots ################
  output$metaplot <- renderPlotly({
      
      req(PlotNameList())
      req(input$meta)
      req(input$cond)
      
      df <- data.frame(
        UMAP_1 = umap[, "harmonyumap_1"],
        UMAP_2 = umap[, "harmonyumap_2"]
      )
      
      df$expr <- as.numeric(expr[input$meta, ])
      df$celltype <- meta$sctype_classification
      df$cluster <- meta$seurat_clusters
      df$treatment <- meta$treatment
      df$experiment <- meta$experiment
      df$orig.ident <- meta$orig.ident
      
      colour_var <- switch(
        input$colour_by,
        "gene" = df$expr,
        "celltype" = df$celltype,
        "cluster" = df$cluster,
        "treatment" = df$treatment,
        "experiment" = df$experiment,
        "orig.ident" = df$orig.ident
      )
      
      df$groupby <- colour_var
# -------------------------------------------------------------------------------
      PlotsList <- list()
      Titles <-list()
      Count <- 1
      
      for (i in PlotNameList()){
        parts <- strsplit(i, ":")[[1]]
        
        column <- tolower(trimws(parts[1]))
        word <- trimws(parts[2])
        
        df_Filter <- df[df[[column]] == word,]
        p <-plot_ly(
            df_Filter,
            x = ~UMAP_1,
            y = ~UMAP_2,
            color = ~groupby,
            type = "scattergl",
            mode = "markers",
            marker = list(size = 3)
          )
        ncol <- 4
# ------ Calculating position of each SubTitle -----------
        tell <- Count
        col <- ((tell - 1) %% ncol) + 1
        row <- ((tell - 1) %/% ncol) + 1
        
        x <- (col - 0.5) / ncol
        y <- 1 - (row - 1) * 0.5 + 0.05
# -----------------------------------------------------        
        PlotsList[[Count]] <- p
        Titles[[Count]] <- list(
          x = x,
          y = y,
          text = word,
          showarrow = FALSE, 
          xanchor = 'center',
          yanchor = 'top',
          font = list(size = 14)
          
        )
        Count <- Count + 1

      }
      
      S <- subplot(PlotsList,
              shareX = TRUE,
              shareY = TRUE,
              nrows = ceiling(Count/4),
              margin = 0.05,
              titleX = TRUE,
              titleY = TRUE
        ) |>layout(annotations= Titles)
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


