# Install the required libraries
library(shiny)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(bslib)
library(thematic)
library(tidyr)
library(plotly)

# Import the data
demo_gene_expression <- read.csv("/home/bg171/srp/final_counts.csv")

title <- tags$img(imageOutput("logo"),
                  'Analysis of Single-cell RNA-seq in Metastatic Breast Cancer', align = "center")

# Define UI theme 
theme <- bs_theme(
  bg = "#FFFFFF",
  fg = "#333333",  # Dark grey for text for readability
  primary = "#e6005c",  # Dark pink for primary buttons and accents 
  base_font = font_google("Lato"),
  heading_font = font_google("Lato")
)

# create the ui
ui <- fluidPage(
  theme = theme,
  
  div(class = "header-container", 
      div(style = "float: left; margin-left: 20px;", 
          h2("Analysis of Single-cell RNA-seq in Metastatic Breast Cancer", style = "center")
      ),
      div(style = "clear: both;")
  ),
  
  # add a link to the git hub account that contains the source code and readME page
  tags$li(class = "dropdown", 
          tags$a(href = "https://github.com/Rouslane-A/RNAseq---Steered-research-project/tree/main", 
                 icon("github"), 
                 "Project Source Code")
  ),
  
  # Create the organisation of the ui
  navlistPanel(
    id = "tabset",
    # project overview tab
    tabPanel("Project Overview", "This project involved re-analysis of the 2023 study by Hamelin et al., titled 'Single-cell Analysis Reveals Inter- and Intratumour Heterogeneity in Metastatic Breast Cancer, which explores transcriptomic changes during breast cancer metastasis. Our approach replicated the methodology used by the authors and simultaneously integrated a recent pipeline for analysing cell-type clustering based on gene expression using the same dataset. Additionally, this web interface was developed to present the outcomes of both the analysis, to facilitate the comparison and interpretation of the findings.",
             imageOutput("methods")),
    # original analysis tab
    tabPanel("Original Analysis",
             fluidRow(
               column(12, 
                      tags$h3("Re-Analysis Results"),
                      selectInput("analysis_selection", "Choose Analysis Results:",
                                  choices = c("Clusters", "Gene Set Enrichment Analysis")), # create a drop down menu for different analyses that were conducted
                      
                      # this panel will only show if 'Clusters' is selected from the drop down mene
                      conditionalPanel(
                        condition = "input.analysis_selection == 'Clusters'",
                        fluidRow(
                          tags$h4("Clustering Results - Initial Clusters"),
                          tags$p("We plotted the initial clustering within the data by the origin and by model."),
                          column(6, imageOutput("F_1")),
                          column(6, imageOutput("F_2")),
                          tags$p("The clustering by model plot (b) shows a good separation of the clusters in five different categories (MDAMB231, PDX1, PDX2, PDX3 and PDX4), with some overlapping between them. On the other hand the clustering by origin (b) 
                                 revealed no clear cluster separation, although the distinct origins are clearly visible."),
                          tags$h4("Clustering Results - Post Analysis"),
                          tags$p("We then carried out our analysis and did some graph-based clustering and cell cycle analysis."),
                          column(6, imageOutput("F_3")),
                          column(6, imageOutput("F_4")),
                          column(6, imageOutput("F_5")),
                          tags$p("25 clusters were found based on the gene expression in cells. 10 clusters were found based on the cell cycle stage cells were in. From the 10 distinct clusters of the cell cycle, there is a predominance of cells in the G1.S/S (orange). The rest of the clusters are not well separated. The 
                                 cluster identity shows a good separation for some clusters and overlapping is also observable between them. 
                                 In the plot, different groups of cells clustered together (X-axis) are compared against the percentage of cells in each cluster (Y-axis). There is a variation in the 
                                 cell cycle stage distribution, with some clusters showing 100% of only one cell cycle or a combination of cycles
                                 (cluster 1, 15, 18, 19, 20-24, 3 and 5-7). There are more clusters in the G2.M, G2.M/M.G1, G2.G2.M and M.G1 cell cycle stages. ")
                        )
                      ),
                      # this panel will only show if 'Gene Set Enrichment Analysis' is chosen from the drop down menu
                      conditionalPanel(
                        condition = "input.analysis_selection == 'Gene Set Enrichment Analysis'",
                        fluidRow(
                          tags$h4("Gene Set Enrichment Analysis"),
                          tags$p("We tried to replicate the researchers Gene Set Enrichment Analysis"),
                          column(6, imageOutput("F6")),
                          column(6, imageOutput("F7")),
                          column(6, imageOutput("F8")),
                          tags$p("Due to many methodolgical difficulties, we managed to produce the heatmaps above. Unfortunately, they are uninterpretable. the first heatmap is cluster 1, the second is cluster 8 and the final heatmap is cluster 24.")
                        )
                      )
               )
             )
             ),
    # Novel analysis tab
    tabPanel("Novel Analysis of Hamelin et al., 2023 study", 
             fluidRow(
               column(12, 
                      tags$h3("Novel Analysis Results"),
                      selectInput("analysis_selection", "Choose Analysis Results:",
                                  choices = c("Clusters", "Gene Set Enrichment Analysis")),
                      # Add your diagrams and text here
                      # Conditional panels to show based on the dropdown selection
                      conditionalPanel(
                        condition = "input.analysis_selection == 'Clusters'",
                        fluidRow(
                          tags$h4("Clustering Results - Initial Clusters"),
                          tags$p("We plotted the initial clustering within the data by the sample origin and by model."),
                          column(6, imageOutput("data_image3")),
                          column(6, imageOutput("data_image4")),
                          tags$p("The clusters by origin showed no specific patterns which reveals the need to analyse 
                                 the cells on the single cell level to understand their gene expression. You can see distinct c
                                 lusters when you cluster by model. This indicates that each model (MDAMB231, PDX1, PDX2, 
                                 PDX3 and PDX4) have distinct gene expression levels as shown by the distinct clusters in the t-SNE graph"),
                          tags$h4("Clustering Results - Post Analysis"),
                          tags$p("We conducted the Seurat Pipeline to analyse the count matrix and carried out graph-based clustering. Below you can
                                 see a t-SNE plot clustering the cells by gene expression and a t-SNE plot clustering the cells by cell cycle stage."),
                          column(6, imageOutput("data_image1")),
                          column(6, imageOutput("data_image2")),
                          tags$p("The clustering plot by gene expression revealed 10 different clousters from the novel analysis. The cell cycle clusters show no distinct
                                 clusters indicating that the cells in each cell cycle stage does not influence the gene expression of the cells."),
                          tags$h4("Cluster gene expression results"),
                          tags$p("Finally, to conclude the clustering analysis, we created a heatmap for each cluster to see if clusters have different gene expression levels for different genes."),
                          column(6, imageOutput("data_image5")),
                          tags$p("The heatmap revealed that every cluster has a high expression for different genes. This can be seen through the pattern of the dot with increased intensity of colour and size.")
                        )
                      ),
                      conditionalPanel(
                        condition = "input.analysis_selection == 'Gene Set Enrichment Analysis'",
                        fluidRow(
                          tags$h4("Gene Set Enrichment Analysis"),
                          tags$p("We first wanted to identify any superclusters that could be present in our data. These superclusters were
                                 then used in a heatmap to see which clusters could be responsioble for the EMT transition or proliferation."),
                          column(6, imageOutput("gsea_supercluster")),
                          column(6, imageOutput("gsea_heatmap")),
                          tags$p("We found 2 superclusters from our data. The heatmap with the EMT-transition and proliferation showed that both superclusters
                                 expressed the proliferation markers significantly highly and supercluster B expressed more EMT-transition markers than supercluster A.")
                        )
                      )
               )
             ),
             # the gene search bar to allow gene name querying
             textInput("gene", "Enter gene name"),
             actionButton("submit", "Submit"),
             plotlyOutput("genePlot")
    ),
    # this is the tab to summarise the project
    tabPanel("Summary of the Project",
             tags$p("More clusters were identified in the Novel Analysis compared to the re-analysis. This could be due to the re-analysis method having difficulties so alternatives were used for the steps that weren't possible to do."),
             tags$p("The problems faced in the re-analysis:"),
             tags$ul(
               tags$li("We could not replicate the study to produce the same count matrix."),
               tags$li("Alternative packages and functions had to be used due to a lack of documentation."),
               tags$li("No information was available on how to drop the duplicates in the researcher's paper."),
               tags$li("We had no knowledge of the annotation file used by the researchers."),
               tags$li("We could not retrieve the SQLite annotation file either."),
               tags$li("We could not replicate the cleaning of the libraries in the count matrix."),
               tags$li("We could not retrieve the genes known to peak in transcription for the cell cycle."),
               tags$li("We were unable to replicate the covariates and regress and get the log-normalized counts using glm.fit.")
             ),
             tags$p("We had to use the count matrix from the research paper to address some of these issues."),
             tags$p("We managed to successfully generate a count matrix but we were unable do the cell cycle, PCA, t-SNE and u-Map analysis on it."),
             tags$p("For the Novel analysis, an alternative pipeline was used (Seurat Pipeline) to the original method."),
             tags$ul(
               tags$li("Fastqc and MultiQC for the quality control"),
               tags$li("TrimGalore! for the filtering and trimming"),
               tags$li("HiSat2 instead of STAR for the indexing and mapping and alignment"),
               tags$li("Seurat pipleine instead of the origional analysis method")
             ),
             tags$p("We managed to generate our own count matrix from the raw sequencing FastQ files"),
             tags$p("We also succesfully completed the Seurat pipeline and the Gene Set Enrichment Analysis.")
    )
  )
)

# Server logic
server <- function(input, output, session) {
  thematic::thematic_shiny()
  
  # Define renderImage function for the images
  
  # t-SNE by gene expression clusters (novel)
  output$data_image1 <- renderImage({
    list(src = "/home/bg171/SRP/www/clustertsne.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by model (novel)
  output$data_image3 <- renderImage({
    list(src = "/home/bg171/SRP/www/tsne by origin.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by cell cycle (novel)
  output$data_image2 <- renderImage({
    list(src = "/home/bg171/SRP/www/cellcycletsne.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by model (novel)
  output$data_image4 <- renderImage({
    list(src = "/home/bg171/SRP/www/tsne by model.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # Dot plot of clusters gene expression (novel)
  output$data_image5 <- renderImage({
    list(src = "/home/bg171/SRP/www/dotplot.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # Heatmap to identify superclusters (novel)
  output$gsea_supercluster <- renderImage({
    list(src = "/home/bg171/SRP/www/superclusterheatmap.png",
         width = "75%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # Heatmap of the GSEA analysis (novel)
  output$gsea_heatmap <- renderImage({
    list(src = "/home/bg171/SRP/www/supercluster_expression_heatmap.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by origin (re-analysis)
  output$F_1 <- renderImage({
    list(src = "/home/bg171/SRP/www/by origin.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by origin (re-analysis)
  output$F_2 <- renderImage({
    list(src = "/home/bg171/SRP/www/bymodelF.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by gene expression clustering (re-analysis)
  output$F_3 <- renderImage({
    list(src = "/home/bg171/SRP/www/graphbased_clustering.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # t-SNE by cell cycle clustering (Re-analysis)
  output$F_4 <- renderImage({
    list(src = "/home/bg171/SRP/www/cc_clustering.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # barplot of cell cycle (re-analysis)
  output$F_5 <- renderImage({
    list(src = "/home/bg171/SRP/www/cc_barplot.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # cluster 1 heatmap (re-analysis)
  output$F6 <- renderImage({
    list(src = "/home/bg171/SRP/www/thumbnail_cluster1.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # cluster 8 heatmap (re-analysis)
  output$F7 <- renderImage({
    list(src = "/home/bg171/SRP/www/cluster8.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # cluster 24 heatmap (re-analysis)
  output$F8 <- renderImage({
    list(src = "/home/bg171/SRP/www/cluster24.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # Methods flowchart comparison
  output$methods <- renderImage({
    list(src = "/home/bg171/SRP/www/method.png",
         width = "100%",
         height = "auto")
  }, deleteFile = FALSE)
  
  # Define server logic for gene plot
  output$genePlot <- renderPlotly({
    req(input$submit)  # Ensure gene submit button is pressed
    
    # plotly graph for the gene query
    gene_data <- demo_gene_expression %>%
      filter(`...1` == input$gene)
    gene_data_long <- gene_data %>%
      gather(key = "Sample", value = "GeneExpression", -"...1")
    print(head(gene_data_long))
    hover_text <- paste("Gene: ", gene_data_long$...1, "<br>",
                        "Sample: ", gene_data_long$Sample, "<br>",
                        "Expression Level: ", gene_data_long$GeneExpression)
    plot_ly(gene_data_long, x = ~"Samples", y = ~GeneExpression, type = 'bar', color = ~Sample,
            text = hover_text, hoverinfo = "text") %>%
      layout(title = paste("Gene Expression of", input$gene, "in Different Samples"),
             xaxis = list(title = "Samples"),
             yaxis = list(title = "Gene Expression Level"))
  })
}

# Run the app
shinyApp(ui, server)
