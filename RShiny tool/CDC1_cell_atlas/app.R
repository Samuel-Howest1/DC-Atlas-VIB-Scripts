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
SeuratObjT <- readRDS("/Users/irc/Desktop/Interschip Bioinformatis 2025-2026/Harmony_Treat_Exp_ADT_Correct_Annot3.rds")
SeuratObjT <- JoinLayers(SeuratObjT,assay = "RNA")
SeuratObjT$cluster_new <- as.character(Idents(SeuratObjT))


# # Define cluster -> cell type mapping
# celltype_map <- c(
#   "0" = "Late Mature",
#   "1" = "Late Immature",
#   "2" = "Early Mature",
#   "3" = "Late Mature",
#   "4" = "Late Mature",
#   "5" = "Late Immature",
#   "6" = "Early Mature",
#   "7" = "Late Mature",
#   "8" = "Early Immature",
#   "9" = "Proliferatingtqble cDC1s",
#   "10" = "Early Immature",
#   "11" = "Proliferating cDC1s",
#   "12" = "Early Immature",
#   "13" = "Proliferating cDC1s",
#   "14" = "Early Immature",
#   "15" = "Proliferating cDC1s",
#   "16" = "Late Immature",
#   "17" = "Late Immature",
#   "18" = "Late Mature",
#   "19" = "Late Mature",
#   "20" = "Early Mature",
#   "21" = "Late Immature",
#   "22" = "Proliferating cDC1s",
#   "23" = "Early Mature",
#   "24" = "Early Mature",
#   "25" = "cDC1s engulfing RBCs",
#   "26" = "Late Immature",
#   "27" = "Late Mature",
#   "28" = "Early Immature",
#   "29" = "Early Immature",
#   "30" = "Late Mature",
#   "31" = "Late Mature",
#   "32" = "Late Mature",
#   "33" = "Early Mature",
#   "34" = "Proliferating cDC1s",
#   "35" = "Late Mature"
# )
# 
# # Add cell type annotation
# seurat_obj$celltype_new <- celltype_map[as.character(seurat_obj$cluster_new)]


Idents(SeuratObjT) <- SeuratObjT$celltype_new
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
metadata <- c("celltype","cluster","treatment","experiment","orig.ident")

genes <- rownames(expr)
########### Avg and Percentage expression of Genes for Dotplot later ##########
# avg <- AverageExpression(SeuratObjT,group.by = "celltype_new",features = genes,layer = "data")$RNA
# 
# pct_group <- meta$celltype_new
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
table_length <- c(seq(1,16))
report_cols <-c(
    "JVE008", "JVE010", "SAM016",
    "SAM05", "SAM06", "SAM2",
    "SAM3", "VBO004", "VBO005",
    "VBO006", "VBO007", "VBO008",
    "VBO009", "VBO010", "VBO011",
    "VBO012"
  )
experiment_map <- c(
  "Test",
  "Test",
  "Final",
  "Final",
  "Notch",
  "LNP_WT",
  "LNP_eLNPs",
  "LNP_pIC_LNPs",
  "LNP_CpG",
  "LNP_pIC",
  "LNP_eLNPs",
  "LNP_pIC_LNPs",
  "LNP_CpG",
  "LNP_pIC",
  "Toxo",
  "Toxo"
)
WT_map <- c(
  "SAM2"   = "WT",
  "SAM3"   = "WT",
  "SAM05"  = "WT",
  "SAM06"  = "WT",
  "SAM016" = "WT",
  "VBO004" = "WT",
  "VBO005" = "Test",
  "VBO006" = "Test",
  "VBO007" = "Test",
  "VBO008" = "Test",
  "VBO009" = "Test",
  "VBO010" = "Test",
  "VBO011" = "Test",
  "VBO012" = "Test",
  "JVE008" = "WT",
  "JVE010" = "Test"
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
      tags$img(src = "lab.png", height = "30px")," Experiment CITEseq:"
    ))
)
tbl <- data.frame(
  Condition = conditions,
  matrix(
    "",
    nrow = length(conditions),
    ncol = length(table_length),
    dimnames = list(NULL, report_cols)
  ),
  check.names = FALSE)
# Fill rows
tbl[1, report_cols] <- report_cols                      # Orig.idents
tbl[2, report_cols] <- WT_map[report_cols]             # Treatment
tbl[3, report_cols] <- experiment_map[report_cols]     # Experiment
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
              overflow: visible;
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
  nav_panel(title = "Gene Plots",
# ------------------ Selectors -----------------------------------
card(   style = "margin-bottom: 5px;",
        max_height = "150px",
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
            
            plotOutput("DimPlotMeta", height = "500")
            
          ),col_widths = 12)
        ),
  
#############################################################################
################# Cell Metdata ######################################################
  nav_panel(title = "Cell Metadata", 
            tags$div( class = "input",
  layout_columns(
                      
            selectizeInput(
              inputId = "meta",
              label = "Select Gene to analysis",
              choices = NULL,
              selected = NULL,
              multiple = F,
              options = list(
              placeholder = 'Type a gene...',
              create = TRUE   # allows typing custom values
              )),
            
            selectizeInput(
              inputId = "cond",
              label = "Select or type the subset",
              choices = condition,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = 'Type a subset...',
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
            ),col_widths= c(4,4,4)),
  
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
          tags$div( class = "input",
            layout_columns(
                    selectizeInput(
                      inputId = "gene_1",
                      label = "Select Gene 1 (X-axis):",
                      choices = NULL,
                      selected = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = 'Type a gene...',
                        create = TRUE   # allows typing custom values
                      )),
                    
                    selectizeInput(
                      inputId = "gene_2",
                      label = "Select Gene 2 (Y-axis):",
                      choices = NULL,
                      selected = NULL,
                      multiple = FALSE,
                      options = list(
                        placeholder = 'Type a gene...',
                        create = TRUE   # allows typing custom values
                      )),
          ),col_widths= c(4,4)),
           
          tags$div(
            actionButton("start_comp", "Generate Comparison"),
            
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
  
  observe({cat("Current tqb is loqdaed")})
  observe({print(input$start)})
  
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
    df$celltype <- meta$celltype_new

    # To male Hover Text of cells
    df$hover <- paste(
      rownames(df),
      "<br>Expression:",
      round(df$expr, 2),
      "<br>Celltype:",df$celltype
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
      celltype = meta$celltype_new
    )
    
    ggplot(df_violin, aes(x = celltype, y= expr,fill = celltype)) +
      geom_violin()+
      scale_fill_brewer(palette = "Set3") +
      theme_bw()+  
      theme(
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
      )
  })
  
  # Info verbeter zodat het duielijkt is 
  # output$DotPlot <- renderPlot({
  #   DotPlot(
  #     SeuratObjT,
  #     features = genes,
  #     group.by = "celltype_new"
  #   )
  
  output$DimPlotMeta <- renderPlot({
    
    req(input$Dimplot)

    df <- data.frame(
      UMAP_1 = umap[, "harmonyumap_1"],
      UMAP_2 = umap[, "harmonyumap_2"]
    )

    df$celltype <- meta$celltype_new
    df$cluster <- meta$cluster_new
    df$treatment <- meta$treatment
    df$experiment <- meta$experiment
    df$orig.ident <- meta$orig.ident

    df$MetaDataColumn <- switch(
      input$Dimplot,
      "celltype" = df$celltype,
      "cluster" = df$cluster,
      "treatment" = df$treatment,
      "experiment" = df$experiment,
      "orig.ident" = df$orig.ident
    )
  ggplot(df,aes(x = UMAP_1, y = UMAP_2, color=MetaDataColumn ))+
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
      df$celltype <- meta$celltype_new
      df$cluster <- meta$cluster_new
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
        
        print(class(p))

        PlotsList[[Count]] <- p

        Count <- Count + 1

      }
      

      subplot(PlotsList,
              shareX = TRUE,
              shareY = TRUE,
              nrows = ceiling(Count/4),
              margin = 0.05
        )
      

    })
  
  
  
  
#################################################################################
Complot <-eventReactive(input$start_comp,{
  input$gene_1
  input$gene_2
})
# ADD COUNT OF HOW MANY CELL PER GENE
  output$Comparison <- renderPlot({
    # FeatureScatter(
    #   SeuratObjT,
    #   feature1 = input$gene_1,
    #   feature2 = input$gene_2
    # )

      df <- data.frame(
        gene1 = as.numeric(expr[input$gene_1,]),
        gene2 = as.numeric(expr[input$gene_2,]),
        color = meta$celltype_new
      )
      
      ggplot(df, aes(x = gene1, y= gene2, color=color))+
        geom_point(size = 0.5)+
        labs(x=input$gene_1, y=input$gene_2 )
      
  })
}
# Run the application 
shinyApp(ui = ui, server = server)


