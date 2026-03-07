# ============================================================================
# BULK RNA-SEQ ANALYSIS SHINY APP - VERSION 2 WITH BATCH CORRECTION
# ============================================================================

# for gsea, uses stats <- sign(comp_df$logFC) * -log10(pvals) * abs(comp_df$logFC) to rank genes/ pathways (better emphasis on biological changes)
# supposedly better batch correction

library(shiny)
library(edgeR)
library(limma)  # For batch correction
library(ggplot2)
library(dplyr)
library(tidyr)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(msigdbr)
library(BiocParallel)
library(pheatmap)
library(grid)
library(gridExtra)
library(ggrepel)
library(writexl)
library(openxlsx)
library(stats)
library(DT)
BiocParallel::register(BiocParallel::SerialParam())

library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(org.Dr.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)

APP_COLORS <- list(
  primary = "#3498db", secondary = "#9b59b6", success = "#27ae60",
  info = "#17a2b8", warning = "#f39c12", danger = "#e74c3c",
  up = "#e74c3c", down = "#3498db", neutral = "#95a5a6"
)

get_orgdb <- function(ref_genome){
  switch(ref_genome,
         "Mus musculus" = org.Mm.eg.db, "Homo sapiens" = org.Hs.eg.db,
         "Rattus norvegicus" = org.Rn.eg.db, "Danio rerio" = org.Dr.eg.db,
         "Drosophila melanogaster" = org.Dm.eg.db, "Caenorhabditis elegans" = org.Ce.eg.db,
         NULL)
}

detect_id_is_ensembl <- function(ids){ mean(grepl("^ENS[A-Z]*G", ids)) > 0.5 }

save_grob_png <- function(grob, file, width=1000, height=800, res=120){
  png(file, width=width, height=height, res=res)
  grid::grid.newpage(); grid::grid.draw(grob); dev.off()
}

# Theme with complete border (for QC, ORA, GSEA, heatmaps, pathway explorer)
theme_publication <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = 13, face = "bold", color = "black"),
      axis.text.y = element_text(size = 13, face = "bold", color = "black"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
}

# Theme for volcano/MA plots (axis lines only, no complete border)
theme_publication_volcano_ma <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = 13, face = "bold", color = "black"),
      axis.text.y = element_text(size = 13, face = "bold", color = "black"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
}

# ============================================================================
# USER INTERFACE
# ============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML(paste0("
      .navbar { background-color: ", APP_COLORS$primary, " !important; }
      .well { background-color: #f8f9fa; border: 1px solid #dee2e6; }
      .btn-primary { background-color: ", APP_COLORS$primary, "; border-color: ", APP_COLORS$primary, "; }
      .btn-primary:hover { background-color: ", APP_COLORS$secondary, "; }
      .btn-success { background-color: ", APP_COLORS$success, "; border-color: ", APP_COLORS$success, "; }
      .btn-info { background-color: ", APP_COLORS$info, "; border-color: ", APP_COLORS$info, "; }
      .plot-container { background-color: white; padding: 20px; border-radius: 5px; 
                       box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 25px; }
      .sidebar-panel { background-color: #f8f9fa; }
      h3, h4 { color: ", APP_COLORS$primary, "; margin-top: 5px; }
      h5 { margin-bottom: 15px; }
      .info-box { background-color: #e7f3ff; border-left: 4px solid ", APP_COLORS$primary, "; 
                 padding: 15px; margin-bottom: 15px; }
      .progress-box { background-color: #d4edda; border-left: 4px solid ", APP_COLORS$success, "; 
                     padding: 12px; margin: 15px 0; font-size: 14px; }
      .main-panel { padding-top: 20px; }
      .batch-panel { background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; 
                    padding: 10px; margin-top: 10px; }
      .batch-header { color: #856404; font-weight: bold; margin-bottom: 10px; }
    ")))
  ),
  
  titlePanel("Bulk RNA-seq Analysis Suite"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar-panel",
      
      wellPanel(
        h4("Upload Data - Batch 1"),
        fileInput("files", "Count Files (TSV)", multiple = TRUE, accept = c(".tsv", ".txt")),
        textInput("group_labels", "Group Labels (comma separated)", 
                  value = "Control,Experimental,Control,Experimental"),
        
        hr(),
        actionButton("add_batch", "Add Batch", class = "btn-info btn-sm", icon = icon("plus")),
        
        # Batch 2 - conditionally shown
        conditionalPanel(
          "input.add_batch > 0",
          div(class = "batch-panel",
              div(class = "batch-header", icon("layer-group"), " Batch 2"),
              fileInput("files_batch2", "Count Files (TSV)", multiple = TRUE, accept = c(".tsv", ".txt")),
              textInput("group_labels_batch2", "Group Labels (comma separated)", value = ""),
              actionButton("add_batch3", "Add Batch 3", class = "btn-info btn-sm", icon = icon("plus"))
          )
        ),
        
        # Batch 3 - conditionally shown
        conditionalPanel(
          "input.add_batch3 > 0",
          div(class = "batch-panel",
              div(class = "batch-header", icon("layer-group"), " Batch 3"),
              fileInput("files_batch3", "Count Files (TSV)", multiple = TRUE, accept = c(".tsv", ".txt")),
              textInput("group_labels_batch3", "Group Labels (comma separated)", value = "")
          )
        ),
        
        hr(),
        
        # Batch correction toggle - only show when multiple batches exist
        conditionalPanel(
          "input.add_batch > 0",
          div(style = "background-color: #e8f4f8; padding: 10px; border-radius: 5px; margin-top: 10px;",
              checkboxInput("apply_batch_correction", 
                            strong("Apply Batch Correction (limma)"), 
                            value = TRUE),
              tags$small("Uses limma::removeBatchEffect on logCPM values")
          )
        ),
        
        hr(),
        selectInput("reference_genome", "Reference Genome", 
                    choices = c("Mus musculus", "Homo sapiens", "Danio rerio",
                                "Drosophila melanogaster", "Rattus norvegicus", "Caenorhabditis elegans"),
                    selected = "Mus musculus")
      ),
      
      wellPanel(
        h4("Gene Filtering"),
        selectInput("gene_filter", "Biotype", 
                    choices = c("All", "Protein Coding", "Long Non-Coding", "MicroRNA"),
                    selected = "Protein Coding"),
        checkboxInput("remove_mito", "Remove Mitochondrial Genes", value = TRUE),
        checkboxInput("remove_rik", "Remove 'Rik' Genes", value = TRUE)
      ),
      
      wellPanel(
        h4("Significance Thresholds"),
        numericInput("fdr_threshold", "FDR Threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("logfc_threshold", "Log2 FC Threshold", value = 1, min = 0, step = 0.1),
        conditionalPanel(
          "input.volcano_pval_type == 'PValue'",
          numericInput("pval_threshold", "P-Value Threshold", value = 0.05, min = 0, max = 1, step = 0.001)
        )
      ),
      
      conditionalPanel(
        "output.data_loaded",
        wellPanel(
          h4("Comparison Selection"),
          selectInput("selected_comparison", "Comparison", choices = NULL)
        )
      ),
      
      wellPanel(
        h4("Gene Set Collection"),
        selectInput("msigdb_collection", NULL,
                    choices = c("Hallmark"="H", "GO Biological Process"="C5_BP", 
                                "GO Cellular Component"="C5_CC", "GO Molecular Function"="C5_MF",
                                "KEGG"="C2_CP:KEGG", "Reactome"="C2_CP:REACTOME"),
                    selected = "H"),
        numericInput("fgsea_topN", "Top N GSEA Pathways", value = 15, min = 1, max = 50, step = 1)
      ),
      
      hr(),
      actionButton("run", "Run Analysis", class = "btn-primary btn-lg btn-block", icon = icon("play")),
      
      conditionalPanel(
        "output.data_loaded",
        hr(),
        downloadButton("download_raw_counts", "Download Raw Counts (.csv)", class = "btn-block"),
        downloadButton("download_results", "Download DE Results (.xlsx)", class = "btn-block"),
        downloadButton("download_siggenes", "Download Sig Genes (.csv)", class = "btn-block"),
        downloadButton("download_gsea_results", "Download GSEA Results (.xlsx)", class = "btn-block")
      )
    ),
    
    mainPanel(
      width = 9,
      class = "main-panel",
      
      uiOutput("progress_box"),
      
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Setup", icon = icon("upload"),
                 div(class = "info-box",
                     h4(icon("info-circle"), " Getting Started"),
                     p("1. Upload your count files (TSV format with gene IDs as rownames)"),
                     p("2. Enter group labels matching file order"),
                     p("3. Optionally add more batches and enable batch correction"),
                     p("4. Configure parameters and click 'Run Analysis'")
                 ),
                 uiOutput("setupSummary"),
                 br(),
                 fluidRow(
                   column(6,
                          div(class = "plot-container",
                              h5("Filtering Summary"),
                              DTOutput("filteringSummary")
                          )
                   ),
                   column(6,
                          div(class = "plot-container",
                              h5("Sample -> Group/Batch Mapping"),
                              DTOutput("sampleGroupMap")
                          )
                   )
                 )
        ),
        
        tabPanel("Quality Control", icon = icon("check-circle"),
                 h3("Data Quality Metrics"),
                 wellPanel(
                   radioButtons("qc_plot_type", "Plot Type:",
                                choices = c("Mitochondrial Fraction" = "mito",
                                            "Library Size" = "boxplot",
                                            "Expression Distribution" = "violin",
                                            "PCA" = "pca"),
                                selected = "pca", inline = TRUE),
                   conditionalPanel(
                     "input.qc_plot_type == 'pca'",
                     checkboxInput("show_ellipses", "Show group ellipses", value = FALSE)
                   )
                 ),
                 div(class = "plot-container",
                     h5(textOutput("qc_plot_title")),
                     plotOutput("qcPlot", height = "600px"),
                     downloadButton("download_qc", "Download", class = "btn-sm")
                 )
        ),
        
        tabPanel("Differential Expression", icon = icon("chart-line"),
                 wellPanel(
                   fluidRow(
                     column(4, radioButtons("de_plot_type", "Plot Type:",
                                            choices = c("Volcano" = "volcano", "MA" = "ma", "Heatmap" = "heatmap"),
                                            selected = "volcano", inline = TRUE)),
                     column(4, 
                            conditionalPanel(
                              "input.de_plot_type == 'volcano'",
                              selectInput("volcano_pval_type", "Y-axis",
                                          choices = c("FDR"="FDR", "P-value"="PValue"), selected = "FDR")
                            ),
                            conditionalPanel(
                              "input.de_plot_type == 'volcano'",
                              numericInput("volcano_ylim_max", "Y-max (optional)", value = NA, min = 0)
                            )
                     ),
                     column(4, 
                            conditionalPanel(
                              "input.de_plot_type == 'heatmap'",
                              radioButtons("sig_heatmap_mode", "Genes:",
                                           choices = c("Union"="union", "Selected"="selected"), 
                                           selected = "union", inline = TRUE)
                            )
                     )
                   )
                 ),
                 div(class = "plot-container",
                     h5(textOutput("de_plot_title")),
                     plotOutput("dePlot", height = "700px"),
                     downloadButton("download_de", "Download", class = "btn-sm")
                 ),
                 fluidRow(
                   column(6,
                          div(class = "plot-container",
                              h5("Significant Gene Counts"),
                              DTOutput("sigCounts")
                          )
                   ),
                   column(6,
                          div(class = "plot-container",
                              h5("Top Results"),
                              DTOutput("resultsTable")
                          )
                   )
                 )
        ),
        
        tabPanel("ORA", icon = icon("list-check"),
                 wellPanel(
                   fluidRow(
                     column(6, radioButtons("ora_plot_type", "Plot Type:", 
                                            choices = c("Dotplot" = "dot", "Barplot" = "bar"),
                                            selected = "dot", inline = TRUE)),
                     column(6, radioButtons("ora_direction", "Direction:", 
                                            choices = c("Both" = "both", "Upregulated" = "up", "Downregulated" = "down"),
                                            selected = "both", inline = TRUE))
                   )
                 ),
                 div(class = "plot-container",
                     h5(textOutput("ora_plot_title")),
                     plotOutput("oraPlot", height = "700px"),
                     downloadButton("download_ora", "Download", class = "btn-sm")
                 )
        ),
        
        tabPanel("GSEA", icon = icon("chart-area"),
                 wellPanel(
                   radioButtons("gsea_plot_type", "Plot Type:", 
                                choices = c("Dotplot" = "dot", "Barplot" = "bar"),
                                selected = "dot", inline = TRUE)
                 ),
                 div(class = "plot-container",
                     h5(textOutput("gsea_plot_title")),
                     plotOutput("gseaPlot", height = "700px"),
                     downloadButton("download_gsea", "Download", class = "btn-sm")
                 )
        ),
        
        tabPanel("Pathway Explorer", icon = icon("microscope"),
                 wellPanel(
                   fluidRow(
                     column(6, selectizeInput("selected_pathways", "Select Pathways", choices = NULL, multiple = TRUE)),
                     column(6, radioButtons("pathway_view", "View:",
                                            choices = c("Comparison Barplot" = "barplot", 
                                                        "Enrichment Traces" = "traces",
                                                        "Gene Heatmap" = "heatmap"),
                                            selected = "barplot", inline = TRUE))
                   ),
                   conditionalPanel(
                     "input.pathway_view == 'barplot'",
                     radioButtons("pathway_barplot_mode", "Comparison Mode:",
                                  choices = c("Selected" = "selected", "vs Control" = "vs_control", 
                                              "All pairwise" = "all_pairwise"),
                                  selected = "selected", inline = TRUE)
                   )
                 ),
                 div(class = "plot-container",
                     h5(textOutput("pathway_plot_title")),
                     plotOutput("pathwayViewPlot", height = "700px"),
                     downloadButton("download_pathway_view", "Download", class = "btn-sm")
                 )
        ),
        
        tabPanel("Custom Genes", icon = icon("star"),
                 wellPanel(
                   fluidRow(
                     column(3, 
                            selectInput("gene_selection_mode", "Gene Selection:",
                                        choices = c("Custom Genes" = "custom", "Top N Genes" = "topn"),
                                        selected = "custom")
                     ),
                     column(3,
                            conditionalPanel(
                              "input.gene_selection_mode == 'custom'",
                              textInput("custom_genes", "Gene Symbols (comma separated)", 
                                        value = "", placeholder = "e.g., Actb, Gapdh, Hprt1")
                            ),
                            conditionalPanel(
                              "input.gene_selection_mode == 'topn'",
                              numericInput("custom_topN", "Number of Genes", value = 20, min = 1, step = 1)
                            )
                     ),
                     column(3,
                            conditionalPanel(
                              "input.gene_selection_mode == 'topn'",
                              radioButtons("custom_topn_mode", "Gene Source:",
                                           choices = c("Selected Comparison" = "selected", "Union (All)" = "union"),
                                           selected = "selected", inline = TRUE)
                            ),
                            conditionalPanel(
                              "input.gene_selection_mode == 'topn' && input.custom_topn_mode == 'selected'",
                              selectInput("custom_topn_comparison", "Comparison:", choices = NULL)
                            )
                     ),
                     column(3,
                            selectInput("sv_metric", "Violin Metric:",
                                        choices = c("logCPM", "per-sample logFC vs reference"),
                                        selected = "logCPM")
                     ),
                     column(3,
                            radioButtons("custom_plot_type", "Plot Type:",
                                         choices = c("Heatmap" = "heatmap", "Dotplot" = "dotplot", "Violin" = "violin"),
                                         selected = "heatmap", inline = TRUE)
                     )
                   )
                 ),
                 div(class = "plot-container",
                     h5(textOutput("custom_plot_title")),
                     plotOutput("customPlot", height = "auto"),
                     downloadButton("download_custom", "Download", class = "btn-sm")
                 )
        )
      )
    )
  )
)


# ============================================================================
# SERVER LOGIC
# ============================================================================

server <- function(input, output, session) {
  
  progress <- NULL
  
  auto_show_rownames <- function(n_genes, threshold = 60) n_genes <= threshold
  
  # Helper function to read count files from a batch
  read_batch_counts <- function(file_paths, file_names) {
    if (is.null(file_paths) || length(file_paths) == 0) return(NULL)
    
    counts_list <- lapply(seq_along(file_paths), function(i) {
      fp <- file_paths[i]
      sample_name <- file_names[i]
      
      df <- tryCatch({
        read.table(fp, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
      }, error = function(e) NULL)
      
      if (is.null(df) || ncol(df) == 0 || any(is.na(colnames(df)))) {
        df2 <- tryCatch({
          tmp <- read.table(fp, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                            col.names = c("gene_id", sample_name), check.names = FALSE)
          rownames(tmp) <- tmp$gene_id
          tmp[, sample_name, drop = FALSE]
        }, error = function(e) NULL)
        if (!is.null(df2)) df <- df2
      }
      
      if (!is.null(df) && ncol(df) == 1 && (is.null(colnames(df)) || is.na(colnames(df)))) {
        colnames(df) <- sample_name
      }
      
      if (is.null(df) || ncol(df) != 1) {
        stop(paste("Could not parse file:", sample_name))
      }
      df
    })
    
    gene_list <- Reduce(intersect, lapply(counts_list, rownames))
    counts_list <- lapply(counts_list, function(df) df[gene_list, , drop = FALSE])
    counts <- do.call(cbind, counts_list)
    colnames(counts) <- file_names
    rownames(counts) <- gene_list
    counts
  }
  
  # Main analysis
  analysisResults <- eventReactive(input$run, {
    progress <<- Progress$new(session, min = 0, max = 10)
    progress$set(message = "Starting analysis...", value = 0)
    on.exit(progress$close())
    
    req(input$files)
    
    progress$set(message = "Reading count files...", value = 1)
    
    # Determine if we'll have multiple batches
    has_batch2 <- input$add_batch > 0 && !is.null(input$files_batch2) && length(input$files_batch2$datapath) > 0
    has_batch3 <- input$add_batch3 > 0 && !is.null(input$files_batch3) && length(input$files_batch3$datapath) > 0
    will_have_multiple_batches <- has_batch2 || has_batch3
    
    # Read Batch 1 - only prefix with B1_ if multiple batches exist
    if (will_have_multiple_batches) {
      batch1_names <- paste0("B1_", input$files$name)
    } else {
      batch1_names <- input$files$name
    }
    counts_batch1 <- read_batch_counts(input$files$datapath, batch1_names)
    group_batch1 <- factor(trimws(unlist(strsplit(input$group_labels, ","))))
    
    validate(need(length(group_batch1) == ncol(counts_batch1),
                  paste("Batch 1: Group labels mismatch:", length(group_batch1), "labels vs", ncol(counts_batch1), "samples")))
    
    batch_info <- rep("Batch1", ncol(counts_batch1))
    all_counts <- counts_batch1
    all_groups <- group_batch1
    
    # Read Batch 2 if exists
    if (has_batch2) {
      batch2_names <- paste0("B2_", input$files_batch2$name)
      counts_batch2 <- read_batch_counts(input$files_batch2$datapath, batch2_names)
      group_batch2 <- factor(trimws(unlist(strsplit(input$group_labels_batch2, ","))))
      
      validate(need(length(group_batch2) == ncol(counts_batch2),
                    paste("Batch 2: Group labels mismatch:", length(group_batch2), "labels vs", ncol(counts_batch2), "samples")))
      
      common_genes <- intersect(rownames(all_counts), rownames(counts_batch2))
      all_counts <- cbind(all_counts[common_genes, ], counts_batch2[common_genes, ])
      all_groups <- factor(c(as.character(all_groups), as.character(group_batch2)))
      batch_info <- c(batch_info, rep("Batch2", ncol(counts_batch2)))
    }
    
    # Read Batch 3 if exists
    if (has_batch3) {
      batch3_names <- paste0("B3_", input$files_batch3$name)
      counts_batch3 <- read_batch_counts(input$files_batch3$datapath, batch3_names)
      group_batch3 <- factor(trimws(unlist(strsplit(input$group_labels_batch3, ","))))
      
      validate(need(length(group_batch3) == ncol(counts_batch3),
                    paste("Batch 3: Group labels mismatch:", length(group_batch3), "labels vs", ncol(counts_batch3), "samples")))
      
      common_genes <- intersect(rownames(all_counts), rownames(counts_batch3))
      all_counts <- cbind(all_counts[common_genes, ], counts_batch3[common_genes, ])
      all_groups <- factor(c(as.character(all_groups), as.character(group_batch3)))
      batch_info <- c(batch_info, rep("Batch3", ncol(counts_batch3)))
    }
    
    counts <- all_counts
    group <- all_groups
    batch <- factor(batch_info)
    has_multiple_batches <- length(unique(batch)) > 1
    
    initial_count <- nrow(counts)
    
    progress$set(message = "Annotating genes...", value = 2)
    
    ref_dataset <- switch(input$reference_genome,
                          "Mus musculus" = "mmusculus_gene_ensembl",
                          "Homo sapiens" = "hsapiens_gene_ensembl",
                          "Danio rerio" = "drerio_gene_ensembl",
                          "Drosophila melanogaster" = "dmelanogaster_gene_ensembl",
                          "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
                          "Caenorhabditis elegans" = "celegans_gene_ensembl")
    
    ensembl <- tryCatch({
      useMart("ensembl", dataset = ref_dataset)
    }, error = function(e) {
      showNotification("Ensembl connection failed - skipping annotation", type = "warning", duration = 3)
      return(NULL)
    })
    
    counts_final <- counts
    after_annotation_count <- nrow(counts)
    
    if (!is.null(ensembl) && (input$gene_filter != "All" || input$remove_mito || input$remove_rik)) {
      attributes <- c('ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'chromosome_name')
      annotations <- tryCatch({
        getBM(attributes = attributes, filters = 'ensembl_gene_id', values = rownames(counts), mart = ensembl)
      }, error = function(e) NULL)
      
      if (!is.null(annotations) && nrow(annotations) > 0) {
        if (input$gene_filter != "All") {
          if (input$gene_filter == "Protein Coding") {
            annotations <- annotations %>% filter(gene_biotype == "protein_coding")
          } else if (input$gene_filter == "Long Non-Coding") {
            annotations <- annotations %>% filter(gene_biotype %in% c("lincRNA", "lncRNA"))
          } else if (input$gene_filter == "MicroRNA") {
            annotations <- annotations %>% filter(gene_biotype == "miRNA")
          }
        }
        
        counts_df <- data.frame(ensembl_gene_id = rownames(counts), counts, check.names = FALSE)
        counts_annotated <- merge(annotations, counts_df, by = 'ensembl_gene_id')
        counts_annotated <- counts_annotated %>% filter(!is.na(external_gene_name) & external_gene_name != "")
        
        if (input$remove_rik) {
          counts_annotated <- counts_annotated %>% filter(!grepl("Rik\\d*$", external_gene_name))
        }
        
        counts_aggregated <- counts_annotated %>%
          group_by(external_gene_name) %>%
          summarise(across(colnames(counts), sum), .groups = "drop")
        
        counts_final <- as.data.frame(counts_aggregated)
        rownames(counts_final) <- counts_final$external_gene_name
        counts_final$external_gene_name <- NULL
      }
    }
    after_annotation_count <- nrow(counts_final)
    
    progress$set(message = "Processing mitochondrial genes...", value = 3)
    
    counts_pre_mito <- counts_final
    mito_syms <- character(0)
    mito_by_chr_syms <- character(0)
    
    if (!is.null(ensembl)) {
      annotations_chr <- tryCatch({
        getBM(attributes = c('ensembl_gene_id','external_gene_name','chromosome_name'),
              filters = 'external_gene_name', values = rownames(counts_pre_mito), mart = ensembl)
      }, error = function(e) NULL)
      
      if (!is.null(annotations_chr) && nrow(annotations_chr) > 0) {
        mito_by_chr_syms <- annotations_chr %>%
          filter(tolower(chromosome_name) %in% c("mt","m","mitochondrion_genome","chrmt","mtdna")) %>%
          pull(external_gene_name)
      }
    }
    
    mito_syms <- union(mito_syms,
                       rownames(counts_pre_mito)[grepl("^MT[-.]", rownames(counts_pre_mito), ignore.case = TRUE)])
    mito_genes <- unique(intersect(union(mito_syms, mito_by_chr_syms), rownames(counts_pre_mito)))
    
    total_counts <- colSums(counts_pre_mito)
    if (length(mito_genes) > 0) {
      mito_counts <- colSums(counts_pre_mito[intersect(rownames(counts_pre_mito), mito_genes), , drop = FALSE])
    } else {
      mito_counts <- setNames(rep(0, ncol(counts_pre_mito)), colnames(counts_pre_mito))
    }
    mito_fraction <- mito_counts / total_counts
    mito_fraction_df <- data.frame(Sample = names(mito_fraction),
                                   MitoFraction = as.numeric(mito_fraction), 
                                   Batch = batch,
                                   row.names = NULL)
    
    if (input$remove_mito && length(mito_genes) > 0) {
      counts_final <- counts_final[!(rownames(counts_final) %in% mito_genes), , drop = FALSE]
    }
    after_mito_count <- nrow(counts_final)
    
    progress$set(message = "Running differential expression...", value = 4)
    
    num_samples <- ncol(counts_final)
    validate(need(length(group) == num_samples,
                  paste("Group labels mismatch:", length(group), "labels vs", num_samples, "samples")))
    
    y <- DGEList(counts = counts_final, group = group)
    y$samples$batch <- batch
    keep <- filterByExpr(y, group = group)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y$samples$lib.size <- colSums(y$counts)
    y <- calcNormFactors(y)
    after_filter_count <- nrow(y$counts)
    
    # Create design matrix - include batch if multiple batches and correction requested
    if (has_multiple_batches && input$apply_batch_correction) {
      design <- model.matrix(~0 + group + batch)
      colnames(design) <- gsub("group", "", colnames(design))
      colnames(design) <- gsub("batch", "", colnames(design))
    } else {
      design <- model.matrix(~0+group)
      colnames(design) <- levels(group)
    }
    
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    
    comparisons <- list()
    grp_levels <- levels(group)
    
    group_cols <- which(colnames(design) %in% grp_levels)
    
    if(length(grp_levels) == 2){
      contrast <- rep(0, ncol(design))
      contrast[group_cols[1]] <- -1
      contrast[group_cols[2]] <- 1
      qlf <- glmQLFTest(fit, contrast = contrast)
      res <- as.data.frame(topTags(qlf, n = nrow(y$counts)))
      res$Gene <- rownames(res)
      comparisons[[ paste0(grp_levels[1], "_vs_", grp_levels[2]) ]] <- res
    } else {
      for(i in 1:(length(grp_levels)-1)){
        for(j in (i+1):length(grp_levels)){
          contrast <- rep(0, ncol(design))
          contrast[group_cols[i]] <- -1
          contrast[group_cols[j]] <- 1
          qlf <- glmQLFTest(fit, contrast = contrast)
          res <- as.data.frame(topTags(qlf, n = nrow(y$counts)))
          res$Gene <- rownames(res)
          comparisons[[ paste0(grp_levels[i], "_vs_", grp_levels[j]) ]] <- res
        }
      }
    }
    
    progress$set(message = "Computing PCA...", value = 5)
    
    logCPM_raw <- cpm(y, log = TRUE)
    
    if (has_multiple_batches && input$apply_batch_correction) {
      logCPM <- removeBatchEffect(logCPM_raw, batch = batch, design = model.matrix(~0 + group))
      
      batch_corrected <- TRUE
    } else {
      logCPM <- logCPM_raw
      batch_corrected <- FALSE
    }
    
    norm_counts_long <- pivot_longer(as.data.frame(logCPM), cols = colnames(logCPM), 
                                     names_to = "Sample", values_to = "Log2_CPM")
    norm_counts_long$Batch <- batch[match(norm_counts_long$Sample, colnames(logCPM))]
    
    pcaData <- prcomp(t(logCPM))
    percentVar <- round(100 * pcaData$sdev^2 / sum(pcaData$sdev^2), 2)
    pca_df <- data.frame(Sample = rownames(pcaData$x), PC1 = pcaData$x[,1],
                         PC2 = pcaData$x[,2], group = y$samples$group, batch = batch)
    
    progress$set(message = "Loading gene sets...", value = 6)
    
    OrgDb <- get_orgdb(input$reference_genome)
    all_genes <- rownames(logCPM)
    use_ensembl_ids <- detect_id_is_ensembl(all_genes)
    
    species_str <- switch(input$reference_genome,
                          "Mus musculus" = "Mus musculus", "Homo sapiens" = "Homo sapiens",
                          "Rattus norvegicus" = "Rattus norvegicus", "Danio rerio" = "Danio rerio",
                          "Drosophila melanogaster" = "Drosophila melanogaster",
                          "Caenorhabditis elegans" = "Caenorhabditis elegans")
    
    collection_choice <- input$msigdb_collection
    msig_res <- NULL
    
    if (startsWith(collection_choice, "C5_")) {
      go_cat <- sub("C5_", "", collection_choice)
      go_subcat <- paste0("GO:", go_cat)
      msig_res <- msigdbr(species = species_str, collection = "C5", subcollection = go_subcat)
    } else if (collection_choice == "H") {
      msig_res <- msigdbr(species = species_str, collection = "H")
    } else if (grepl(":", collection_choice)) {
      parts <- strsplit(collection_choice, "_")[[1]]
      main_coll <- parts[1]
      sub_parts <- strsplit(parts[2], ":")[[1]]
      sub_coll <- sub_parts[1]
      sub_cat <- if(length(sub_parts) > 1) sub_parts[2] else NULL
      
      if (!is.null(sub_cat)) {
        msig_res <- msigdbr(species = species_str, collection = main_coll,
                            subcollection = paste0(sub_coll, ":", sub_cat))
      } else {
        msig_res <- msigdbr(species = species_str, collection = main_coll, subcollection = sub_coll)
      }
    } else {
      msig_res <- msigdbr(species = species_str, collection = collection_choice)
    }
    
    pathways_list <- list()
    if (!is.null(msig_res) && nrow(msig_res) > 0) {
      if (use_ensembl_ids && !is.null(OrgDb)) {
        sym_keys <- unique(msig_res$gene_symbol)
        sym2ens <- AnnotationDbi::mapIds(OrgDb, keys = sym_keys, keytype = "SYMBOL",
                                         column = "ENSEMBL", multiVals = "first")
        msig_res$mapped <- unname(sym2ens[msig_res$gene_symbol])
        pathways_list <- split(msig_res$mapped, msig_res$gs_name)
      } else {
        pathways_list <- split(msig_res$gene_symbol, msig_res$gs_name)
      }
      pathways_list <- lapply(pathways_list, function(v) unique(na.omit(v)))
      pathways_list <- pathways_list[vapply(pathways_list, length, 1L) > 0]
    }
    
    enrich_list <- list()
    fgsea_list <- list()
    
    progress$set(message = "Running pathway enrichment...", value = 7)
    
    # Build TERM2GENE table for enricher() - works for any MSigDB collection
    term2gene <- NULL
    if (!is.null(msig_res) && nrow(msig_res) > 0) {
      if (use_ensembl_ids) {
        term2gene <- msig_res %>% dplyr::select(gs_name, mapped) %>% dplyr::filter(!is.na(mapped)) %>%
          dplyr::rename(term = gs_name, gene = mapped)
      } else {
        term2gene <- msig_res %>% dplyr::select(gs_name, gene_symbol) %>%
          dplyr::rename(term = gs_name, gene = gene_symbol)
      }
    }
    
    for(comp_name in names(comparisons)){
      comp_df <- comparisons[[comp_name]]
      
      sig_up <- comp_df %>% filter(FDR <= input$fdr_threshold & logFC >= input$logfc_threshold) %>% pull(Gene)
      sig_down <- comp_df %>% filter(FDR <= input$fdr_threshold & logFC <= -input$logfc_threshold) %>% pull(Gene)
      
      enrich_up <- enrich_down <- NULL
      
      if (!is.null(term2gene) && nrow(term2gene) > 0) {
        universe_genes <- if (use_ensembl_ids && !is.null(msig_res)) {
          unique(na.omit(msig_res$mapped))
        } else {
          all_genes
        }
        
        if (length(sig_up) > 0) {
          enrich_up <- tryCatch({
            enricher(gene = sig_up, TERM2GENE = term2gene, universe = universe_genes,
                     pAdjustMethod = "BH", qvalueCutoff = 0.05)
          }, error = function(e) NULL)
        }
        
        if (length(sig_down) > 0) {
          enrich_down <- tryCatch({
            enricher(gene = sig_down, TERM2GENE = term2gene, universe = universe_genes,
                     pAdjustMethod = "BH", qvalueCutoff = 0.05)
          }, error = function(e) NULL)
        }
      }
      
      enrich_list[[comp_name]] <- list(up = enrich_up, down = enrich_down)
      
      fgseaRes <- NULL
      if (length(pathways_list) > 0) {
        pvals <- comp_df$PValue
        if (all(is.na(pvals))) {
          stats <- rep(NA_real_, length(pvals))
        } else {
          min_nonzero <- min(pvals[pvals > 0 & is.finite(pvals)], na.rm = TRUE)
          pvals[!is.finite(pvals) | pvals <= 0] <- min_nonzero
          stats <- sign(comp_df$logFC) * -log10(pvals) * abs(comp_df$logFC)
        }
        names(stats) <- comp_df$Gene
        stats <- stats[is.finite(stats)]
        stats <- sort(stats, decreasing = TRUE)
        
        if (length(stats) > 1) {
          fgseaRes <- tryCatch({
            fgseaMultilevel(pathways = pathways_list, stats = stats)
          }, error = function(e) NULL)
          if(!is.null(fgseaRes)){
            fgseaRes <- fgseaRes %>% arrange(padj)
          }
        }
      }
      fgsea_list[[comp_name]] <- fgseaRes
    }
    
    progress$set(message = "Finalizing...", value = 9)
    
    sig_counts <- data.frame(
      Comparison = names(comparisons),
      Upregulated = sapply(comparisons, function(df) 
        sum(df$FDR <= input$fdr_threshold & df$logFC >= input$logfc_threshold)),
      Downregulated = sapply(comparisons, function(df) 
        sum(df$FDR <= input$fdr_threshold & df$logFC <= -input$logfc_threshold))
    )
    
    sig_genes_by_comp <- lapply(comparisons, function(df) {
      df %>% filter(FDR <= input$fdr_threshold & abs(logFC) >= input$logfc_threshold) %>%
        pull(Gene) %>% unique()
    })
    all_sig_genes <- unique(unlist(sig_genes_by_comp))
    
    filteringSummary <- data.frame(
      Step = c("Initial", "After Annotation", "After Mito", "After FilterByExpr"),
      GeneCount = c(initial_count, after_annotation_count, after_mito_count, after_filter_count)
    )
    
    progress$set(message = "Complete!", value = 10)
    showNotification("Analysis complete!", type = "message", duration = 5)
    
    list(counts_final = counts_final, counts_raw = counts_pre_mito,
         mito_fraction_df = mito_fraction_df, y = y,
         comparisons = comparisons, logCPM = logCPM, logCPM_raw = logCPM_raw,
         norm_counts_long = norm_counts_long,
         pca_df = pca_df, percentVar = percentVar, enrich_list = enrich_list,
         fgsea_list = fgsea_list, pathways_list = pathways_list,
         sig_genes_union = all_sig_genes, sig_genes_by_comp = sig_genes_by_comp,
         filteringSummary = filteringSummary, sig_counts = sig_counts,
         batch = batch, has_multiple_batches = has_multiple_batches,
         batch_corrected = batch_corrected)
  })
  
  observeEvent(analysisResults(), {
    comp_names <- names(analysisResults()$comparisons)
    updateSelectInput(session, "selected_comparison", choices = comp_names, selected = comp_names[1])
    updateSelectInput(session, "custom_topn_comparison", choices = comp_names, selected = comp_names[1])
  })
  
  observeEvent(input$selected_comparison, {
    req(analysisResults())
    fgseaRes <- analysisResults()$fgsea_list[[input$selected_comparison]]
    if (!is.null(fgseaRes) && nrow(fgseaRes) > 0) {
      pathway_choices <- unique(fgseaRes$pathway)
      updateSelectizeInput(session, "selected_pathways", choices = pathway_choices, 
                           selected = head(pathway_choices, 3), server = TRUE)
    } else {
      updateSelectizeInput(session, "selected_pathways", choices = NULL, selected = NULL, server = TRUE)
    }
  })
  
  selectedDE <- reactive({
    req(analysisResults(), input$selected_comparison)
    res <- analysisResults()$comparisons[[input$selected_comparison]]
    if(!("Gene" %in% colnames(res))) res$Gene <- rownames(res)
    res
  })
  
  selectedGenes <- reactive({
    req(analysisResults())
    mode <- input$gene_selection_mode
    
    if (mode == "custom") {
      if (is.null(input$custom_genes) || !nzchar(trimws(input$custom_genes))) return(character(0))
      input_genes <- unique(trimws(unlist(strsplit(input$custom_genes, ","))))
      all_genes <- rownames(analysisResults()$logCPM)
      matched_genes <- all_genes[toupper(all_genes) %in% toupper(input_genes)]
      return(matched_genes)
    } else if (mode == "topn") {
      req(input$custom_topN)
      n <- as.integer(input$custom_topN)
      validate(need(!is.na(n) && n > 0, "Please set Top N genes > 0."))
      comps <- analysisResults()$comparisons
      if (input$custom_topn_mode == "union") {
        merged <- dplyr::bind_rows(lapply(comps, function(df) {
          if (!("Gene" %in% colnames(df))) df$Gene <- rownames(df)
          df
        }))
        merged_rank <- merged %>%
          dplyr::group_by(Gene) %>%
          dplyr::summarise(minFDR = min(FDR, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(minFDR)
        return(head(merged_rank$Gene, n))
      } else {
        req(input$custom_topn_comparison)
        comp_df <- comps[[input$custom_topn_comparison]]
        if(!("Gene" %in% colnames(comp_df))) comp_df$Gene <- rownames(comp_df)
        return(head(comp_df %>% dplyr::arrange(FDR) %>% dplyr::pull(Gene), n))
      }
    }
    return(character(0))
  })
  
  comparisonGroups <- reactive({
    req(input$selected_comparison)
    parts <- strsplit(input$selected_comparison, "_vs_")[[1]]
    if (length(parts) != 2) return(list(ref = NA_character_, alt = NA_character_))
    list(ref = parts[1], alt = parts[2])
  })
  
  output$data_loaded <- reactive({ !is.null(analysisResults()) })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  output$progress_box <- renderUI({
    req(analysisResults())
    batch_text <- if (analysisResults()$has_multiple_batches) {
      if (analysisResults()$batch_corrected) {
        paste0(" | ", length(unique(analysisResults()$batch)), " batches (corrected)")
      } else {
        paste0(" | ", length(unique(analysisResults()$batch)), " batches (uncorrected)")
      }
    } else ""
    
    div(class = "progress-box",
        strong(icon("check-circle"), " Analysis Complete"),
        " - ", ncol(analysisResults()$y$counts), " samples, ",
        nrow(analysisResults()$y$counts), " genes, ",
        length(analysisResults()$comparisons), " comparisons",
        batch_text
    )
  })
  
  output$setupSummary <- renderUI({
    req(analysisResults())
    batch_info <- if (analysisResults()$has_multiple_batches) {
      p(strong("Batches:"), length(unique(analysisResults()$batch)), 
        if (analysisResults()$batch_corrected) " (batch corrected)" else " (no correction)")
    } else NULL
    
    div(class = "plot-container",
        h4("Analysis Summary"),
        p(strong("Samples:"), ncol(analysisResults()$counts_final)),
        p(strong("Groups:"), paste(unique(analysisResults()$y$samples$group), collapse = ", ")),
        batch_info,
        p(strong("Genes retained:"), nrow(analysisResults()$y$counts)),
        p(strong("Comparisons:"), length(analysisResults()$comparisons))
    )
  })
  
  output$filteringSummary <- renderDT({
    req(analysisResults())
    datatable(analysisResults()$filteringSummary, options = list(dom = 't'), rownames = FALSE)
  })
  
  output$sampleGroupMap <- renderDT({
    req(analysisResults())
    if (analysisResults()$has_multiple_batches) {
      df <- data.frame(Sample = colnames(analysisResults()$counts_final),
                       Group = as.character(analysisResults()$y$samples$group),
                       Batch = as.character(analysisResults()$batch))
    } else {
      df <- data.frame(Sample = colnames(analysisResults()$counts_final),
                       Group = as.character(analysisResults()$y$samples$group))
    }
    datatable(df, options = list(dom = 't', pageLength = 20), rownames = FALSE)
  })
  
  output$resultsTable <- renderDT({
    req(selectedDE())
    datatable(head(selectedDE(), 20), options = list(dom = 't', scrollX = TRUE), rownames = FALSE)
  })
  
  output$sigCounts <- renderDT({
    req(analysisResults())
    datatable(analysisResults()$sig_counts, options = list(dom = 't'), rownames = FALSE)
  })
  
  # QC plot title
  output$qc_plot_title <- renderText({
    switch(input$qc_plot_type,
           "mito" = "Mitochondrial Fraction",
           "boxplot" = "Library Size Distribution",
           "violin" = "Expression Distribution",
           "pca" = "PCA Plot")
  })
  
  # QC plot
  output$qcPlot <- renderPlot({
    req(analysisResults())
    
    if (input$qc_plot_type == "mito") {
      mito_df <- analysisResults()$mito_fraction_df
      if (analysisResults()$has_multiple_batches) {
        p <- ggplot(mito_df, aes(x = Sample, y = MitoFraction, fill = Batch)) +
          geom_bar(stat = "identity", alpha = 0.8) +
          scale_fill_brewer(palette = "Set2")
      } else {
        p <- ggplot(mito_df, aes(x = Sample, y = MitoFraction)) +
          geom_bar(stat = "identity", alpha = 0.8, fill = APP_COLORS$primary)
      }
      p <- p +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = APP_COLORS$danger) +
        ylim(0, 1) +
        theme_publication() +
        labs(x = "Sample", y = "Mitochondrial Fraction") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      
    } else if (input$qc_plot_type == "boxplot") {
      df <- analysisResults()$norm_counts_long
      if (analysisResults()$has_multiple_batches) {
        p <- ggplot(df, aes(x = Sample, y = Log2_CPM, fill = Batch)) +
          geom_boxplot(alpha = 0.6, outlier.color = APP_COLORS$danger) +
          scale_fill_brewer(palette = "Set2")
      } else {
        p <- ggplot(df, aes(x = Sample, y = Log2_CPM)) +
          geom_boxplot(alpha = 0.6, outlier.color = APP_COLORS$danger, fill = APP_COLORS$primary)
      }
      p <- p +
        theme_publication() + labs(x = "Sample", y = "Log2 CPM") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      
    } else if (input$qc_plot_type == "violin") {
      df <- analysisResults()$norm_counts_long
      if (analysisResults()$has_multiple_batches) {
        p <- ggplot(df, aes(x = Sample, y = Log2_CPM, fill = Batch)) +
          geom_violin(alpha = 0.6) +
          geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
          scale_fill_brewer(palette = "Set2")
      } else {
        p <- ggplot(df, aes(x = Sample, y = Log2_CPM)) +
          geom_violin(alpha = 0.6, fill = APP_COLORS$primary) +
          geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA)
      }
      p <- p +
        theme_publication() + labs(x = "Sample", y = "Log2 CPM") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      
    } else if (input$qc_plot_type == "pca") {
      df <- analysisResults()$pca_df
      percentVar <- analysisResults()$percentVar
      has_batches <- analysisResults()$has_multiple_batches
      
      if (has_batches) {
        # Color by group, shape by batch
        p <- ggplot(df, aes(x = PC1, y = PC2, color = group, shape = batch, label = Sample)) +
          geom_point(size = 5, alpha = 0.8) + 
          ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
          scale_color_brewer(palette = "Set2") +
          scale_shape_manual(values = c(16, 17, 15, 18, 8, 3)) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          theme_publication() +
          theme(legend.title = element_blank(), legend.position = "top") +
          guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
      } else {
        # Just color by group
        p <- ggplot(df, aes(x = PC1, y = PC2, color = group, label = Sample)) +
          geom_point(size = 5, alpha = 0.8) + 
          ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
          scale_color_brewer(palette = "Set2") +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          theme_publication() +
          theme(legend.title = element_blank(), legend.position = "top")
      }
      
      if (input$show_ellipses) {
        p <- p + stat_ellipse(aes(color = group), type = "t", level = 0.95, linewidth = 1)
      }
      print(p)
    }
  })
  
  # DE plot title
  output$de_plot_title <- renderText({
    switch(input$de_plot_type,
           "volcano" = "Volcano Plot",
           "ma" = "MA Plot",
           "heatmap" = "Significant Genes Heatmap")
  })
  
  # DE plot
  output$dePlot <- renderPlot({
    req(analysisResults())
    
    if (input$de_plot_type == "volcano") {
      tryCatch({
        df <- selectedDE()
        pcol <- input$volcano_pval_type
        pthresh <- if (pcol=="FDR") input$fdr_threshold else input$pval_threshold
        logfc_thresh <- max(0.01, input$logfc_threshold)
        
        df$Significance <- ifelse(df[[pcol]] <= pthresh & abs(df$logFC) >= logfc_thresh,
                                  "Significant", "Not Significant")
        df$yval <- -log10(df[[pcol]])
        df <- df[is.finite(df$yval) & is.finite(df$logFC), ]
        sig_df <- df[df$Significance == "Significant", ]
        
        labels_df <- data.frame()
        if (nrow(sig_df) > 60) {
          left_labels <- head(sig_df %>% filter(logFC < 0) %>% arrange(!!rlang::sym(pcol), desc(abs(logFC))), 30)
          right_labels <- head(sig_df %>% filter(logFC > 0) %>% arrange(!!rlang::sym(pcol), desc(abs(logFC))), 30)
          labels_df <- bind_rows(left_labels, right_labels)
        } else {
          labels_df <- sig_df
        }
        
        p <- ggplot(df, aes(x = logFC, y = yval, color = Significance)) +
          geom_point(alpha = 0.6, size = 2) +
          scale_color_manual(values = c("Significant" = APP_COLORS$danger, "Not Significant" = "black"))
        
        if (nrow(labels_df) > 0) {
          p <- p + ggrepel::geom_text_repel(data = labels_df, aes(label = Gene), 
                                            size = 3, max.overlaps = 100, color = "black")
        }
        
        p <- p +
          geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "black") +
          {if(logfc_thresh > 0) geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") else NULL} +
          coord_cartesian(ylim = if (is.na(input$volcano_ylim_max)) NULL else c(0, input$volcano_ylim_max)) +
          labs(title = input$selected_comparison, x = "Log2 Fold Change", y = paste0("-Log10(", pcol, ")")) +
          theme_publication_volcano_ma() + theme(legend.position = "top")
        print(p)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error generating volcano plot:\n", e$message), cex = 1.2, col = "red")
      })
      
    } else if (input$de_plot_type == "ma") {
      df <- selectedDE()
      logfc_thresh <- max(0.01, input$logfc_threshold)
      df$Significance <- ifelse(df$FDR <= input$fdr_threshold & abs(df$logFC) >= logfc_thresh,
                                "Significant", "Not Significant")
      avg_logCPM <- rowMeans(analysisResults()$logCPM)
      df$avg_logCPM <- avg_logCPM[match(df$Gene, rownames(analysisResults()$logCPM))]
      
      ggplot(df, aes(x = avg_logCPM, y = logFC, color = Significance)) +
        geom_point(alpha = 0.6, size = 2) +
        scale_color_manual(values = c("Significant" = APP_COLORS$danger, "Not Significant" = "black")) +
        {if(logfc_thresh > 0) geom_hline(yintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") else NULL} +
        theme_publication_volcano_ma() + 
        labs(title = input$selected_comparison, x = "Average Log2 CPM", y = "Log2 Fold Change") +
        theme(legend.position = "top")
      
    } else if (input$de_plot_type == "heatmap") {
      logCPM <- analysisResults()$logCPM
      mode <- input$sig_heatmap_mode
      sig_genes <- if (identical(mode, "union")) {
        analysisResults()$sig_genes_union
      } else {
        req(input$selected_comparison)
        analysisResults()$sig_genes_by_comp[[input$selected_comparison]]
      }
      sig_logCPM <- logCPM[rownames(logCPM) %in% sig_genes, , drop = FALSE]
      if(nrow(sig_logCPM) > 1){
        n_genes <- nrow(sig_logCPM)
        fontsize_row <- max(6, 10 - (n_genes/20))
        if (analysisResults()$has_multiple_batches) {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group,
                                       Batch = analysisResults()$batch)
        } else {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
        }
        rownames(annotation_col) <- colnames(sig_logCPM)
        main_title <- if (identical(mode, "union")) "Significant Genes (All)" else paste("Sig Genes:", input$selected_comparison)
        pheatmap(as.matrix(sig_logCPM), annotation_col = annotation_col, scale = "row",
                 main = main_title, fontsize_row = fontsize_row, fontsize_col = 11,
                 show_rownames = auto_show_rownames(nrow(sig_logCPM)),
                 color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
      } else {
        plot.new(); text(0.5, 0.5, "Not enough significant genes", cex = 1.2)
      }
    }
  })
  
  # ORA plot title
  output$ora_plot_title <- renderText({
    if (input$ora_plot_type == "dot") "ORA Dotplot" else "ORA Barplot"
  })
  
  # ORA plot
  output$oraPlot <- renderPlot({
    req(analysisResults())
    enrich_results <- analysisResults()$enrich_list[[input$selected_comparison]]
    
    plots <- list()
    show_up <- input$ora_direction %in% c("both", "up")
    show_down <- input$ora_direction %in% c("both", "down")
    
    if(show_up && !is.null(enrich_results$up) && nrow(as.data.frame(enrich_results$up)) > 0){
      if (input$ora_plot_type == "dot") {
        p_up <- dotplot(enrich_results$up, showCategory = 20) + 
          ggtitle(paste("Upregulated:", input$selected_comparison)) +
          scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + 
          theme_publication()
      } else {
        p_up <- barplot(enrich_results$up, showCategory = 20) + 
          ggtitle(paste("Upregulated:", input$selected_comparison)) +
          scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + 
          theme_publication()
      }
      plots[["up"]] <- p_up
    }
    
    if(show_down && !is.null(enrich_results$down) && nrow(as.data.frame(enrich_results$down)) > 0){
      if (input$ora_plot_type == "dot") {
        p_down <- dotplot(enrich_results$down, showCategory = 20) + 
          ggtitle(paste("Downregulated:", input$selected_comparison)) +
          scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + 
          theme_publication()
      } else {
        p_down <- barplot(enrich_results$down, showCategory = 20) + 
          ggtitle(paste("Downregulated:", input$selected_comparison)) +
          scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + 
          theme_publication()
      }
      plots[["down"]] <- p_down
    }
    
    if(length(plots) == 0){
      plot.new(); text(0.5, 0.5, "No enriched terms found", cex = 1.2)
    } else if(length(plots) == 1){
      print(plots[[1]])
    } else {
      grid.arrange(plots[["up"]], plots[["down"]], ncol = 1, heights = c(1, 1))
    }
  }, height = function() {
    if (input$ora_direction == "both") return(1300) else return(700)
  })
  
  # GSEA plot title
  output$gsea_plot_title <- renderText({
    if (input$gsea_plot_type == "dot") "GSEA Dotplot" else "GSEA Barplot"
  })
  
  # GSEA plot
  output$gseaPlot <- renderPlot({
    req(analysisResults(), input$fgsea_topN)
    fgseaRes <- analysisResults()$fgsea_list[[input$selected_comparison]]
    if (is.null(fgseaRes) || nrow(fgseaRes) == 0) {
      plot.new(); text(0.5, 0.5, "No significant GSEA pathways", cex = 1.2)
      return()
    }
    topN <- max(1, as.integer(input$fgsea_topN))
    fg_up <- fgseaRes %>% filter(padj < 0.05, NES > 0) %>% arrange(padj) %>% head(topN)
    fg_down <- fgseaRes %>% filter(padj < 0.05, NES < 0) %>% arrange(padj) %>% head(topN)
    topPathways <- bind_rows(fg_up, fg_down)
    if (nrow(topPathways) == 0) {
      plot.new(); text(0.5, 0.5, "No significant GSEA pathways", cex = 1.2)
      return()
    }
    
    if (input$gsea_plot_type == "dot") {
      ggplot(topPathways, aes(x = NES, y = reorder(pathway, NES), size = size, color = -log10(padj))) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) +
        theme_publication() + 
        labs(title = input$selected_comparison, x = "NES", y = "Pathway") +
        theme(axis.text.y = element_text(size = 11, face = "bold"))
    } else {
      ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
        geom_col(alpha = 0.8) + coord_flip() +
        labs(title = input$selected_comparison, x = "Pathway", y = "NES") +
        scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) +
        theme_publication() +
        theme(axis.text.y = element_text(size = 11, face = "bold"))
    }
  }, height = 900)
  
  # Pathway Explorer plot
  output$pathway_plot_title <- renderText({
    switch(input$pathway_view,
           "barplot" = "Pathway Enrichment Comparison",
           "traces" = "Enrichment Traces",
           "heatmap" = "Pathway Genes Heatmap")
  })
  
  output$pathwayViewPlot <- renderPlot({
    req(analysisResults(), input$selected_pathways)
    
    if (input$pathway_view == "barplot") {
      mode <- input$pathway_barplot_mode
      fgsea_list <- analysisResults()$fgsea_list
      comp_names <- names(fgsea_list)
      if (mode == "selected") {
        req(input$selected_comparison)
        comp_names <- input$selected_comparison
      } else if (mode == "vs_control") {
        grp_levels <- levels(analysisResults()$y$samples$group)
        control_group <- grp_levels[1]
        comp_names <- grep(paste0("^", control_group, "_vs_|_vs_", control_group, "$"), comp_names, value = TRUE)
      }
      pathway_data <- lapply(comp_names, function(comp) {
        fgsea_res <- fgsea_list[[comp]]
        if (is.null(fgsea_res)) return(NULL)
        parts <- strsplit(comp, "_vs_")[[1]]
        fgsea_res %>% filter(pathway %in% input$selected_pathways) %>%
          mutate(Comparison = comp, Reference = parts[1], Group = parts[2],
                 sig_label = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))
      })
      pathway_data <- bind_rows(pathway_data)
      if (nrow(pathway_data) == 0) {
        plot.new(); text(0.5, 0.5, "No pathway data", cex = 1.2)
        return()
      }
      p <- ggplot(pathway_data, aes(x = pathway, y = NES, fill = Group)) +
        geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
        geom_col(position = position_dodge(width = 0.8), color = NA, width = 0.75, alpha = 0.8) +
        geom_text(aes(label = sig_label, y = NES + sign(NES) * max(abs(NES)) * 0.08),
                  position = position_dodge(width = 0.8), size = 5, fontface = "bold") +
        coord_flip() + scale_fill_brewer(palette = "Set2") +
        theme_publication() + 
        labs(x = "Pathway", y = "Normalized Enrichment Score (NES)") +
        theme(legend.position = "right")
      if (mode == "all_pairwise") {
        p <- p + facet_wrap(~ paste0("Ref: ", Reference), ncol = 2)
      }
      print(p)
      
    } else if (input$pathway_view == "traces") {
      comp_df <- analysisResults()$comparisons[[input$selected_comparison]]
      pvals <- comp_df$PValue
      min_nonzero <- min(pvals[pvals > 0 & is.finite(pvals)], na.rm = TRUE)
      pvals[!is.finite(pvals) | pvals <= 0] <- min_nonzero
      stats <- sign(comp_df$logFC) * -log10(pvals) * abs(comp_df$logFC)
      names(stats) <- comp_df$Gene
      stats <- stats[is.finite(stats)]
      stats <- sort(stats, decreasing = TRUE)
      plots <- lapply(input$selected_pathways, function(pw) {
        if (pw %in% names(analysisResults()$pathways_list)) {
          plotEnrichment(analysisResults()$pathways_list[[pw]], stats) +
            labs(title = pw) + theme_publication()
        } else NULL
      })
      plots <- Filter(Negate(is.null), plots)
      if(length(plots) == 1){
        plots[[1]]
      } else if(length(plots) > 1){
        do.call(grid.arrange, c(plots, ncol = 1))
      } else {
        plot.new(); text(0.5, 0.5, "No pathway selected", cex = 1.2)
      }
      
    } else if (input$pathway_view == "heatmap") {
      logCPM <- analysisResults()$logCPM
      pathways_list <- analysisResults()$pathways_list
      
      gene_to_pathway <- list()
      for (pw in input$selected_pathways) {
        if (pw %in% names(pathways_list)) {
          genes <- pathways_list[[pw]]
          for (g in genes) {
            if (g %in% rownames(logCPM)) {
              if (is.null(gene_to_pathway[[g]])) {
                gene_to_pathway[[g]] <- pw
              } else {
                gene_to_pathway[[g]] <- paste(gene_to_pathway[[g]], pw, sep = "; ")
              }
            }
          }
        }
      }
      
      selected_genes <- names(gene_to_pathway)
      selected_logCPM <- logCPM[selected_genes, , drop = FALSE]
      
      if(nrow(selected_logCPM) > 1){
        if (analysisResults()$has_multiple_batches) {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group,
                                       Batch = analysisResults()$batch)
        } else {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
        }
        rownames(annotation_col) <- colnames(selected_logCPM)
        
        annotation_row <- data.frame(Pathway = unlist(gene_to_pathway[rownames(selected_logCPM)]))
        rownames(annotation_row) <- rownames(selected_logCPM)
        
        pheatmap(as.matrix(selected_logCPM), 
                 annotation_col = annotation_col,
                 annotation_row = annotation_row,
                 scale = "row",
                 fontsize_row = 9, fontsize_col = 11,
                 show_rownames = auto_show_rownames(nrow(selected_logCPM)),
                 color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
      } else {
        plot.new(); text(0.5, 0.5, "Not enough genes", cex = 1.2)
      }
    }
  })
  
  # Custom genes plot
  output$custom_plot_title <- renderText({
    switch(input$custom_plot_type,
           "heatmap" = "Gene Heatmap",
           "dotplot" = "Gene Expression Dotplot",
           "violin" = "Gene Expression by Group")
  })
  
  output$customPlot <- renderPlot({
    req(analysisResults())
    genes <- selectedGenes()
    
    if (length(genes) == 0) {
      plot.new()
      text(0.5, 0.5, "Select genes using the options above", cex = 1.2)
      return()
    }
    
    logCPM <- analysisResults()$logCPM
    
    if (input$custom_plot_type == "heatmap") {
      selected_logCPM <- logCPM[intersect(rownames(logCPM), genes), , drop = FALSE]
      if(nrow(selected_logCPM) > 1){
        if (analysisResults()$has_multiple_batches) {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group,
                                       Batch = analysisResults()$batch)
        } else {
          annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
        }
        rownames(annotation_col) <- colnames(selected_logCPM)
        pheatmap(as.matrix(selected_logCPM), annotation_col = annotation_col, scale = "row",
                 fontsize_row = 10, fontsize_col = 11,
                 show_rownames = auto_show_rownames(nrow(selected_logCPM)),
                 color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
      } else if(nrow(selected_logCPM) == 1){
        plot.new(); text(0.5, 0.5, "Only one gene\nHeatmap requires >= 2 genes", cex = 1.2)
      } else {
        plot.new(); text(0.5, 0.5, "No matching genes found", cex = 1.2)
      }
      
    } else if (input$custom_plot_type == "dotplot") {
      dot_logCPM <- logCPM[intersect(rownames(logCPM), genes), , drop = FALSE]
      if(nrow(dot_logCPM) > 0){
        dot_df <- as.data.frame(dot_logCPM)
        dot_df$Gene <- rownames(dot_df)
        dot_df_long <- pivot_longer(dot_df, cols = -Gene, names_to = "Sample", values_to = "Expression")
        ggplot(dot_df_long, aes(x = Sample, y = Gene, size = Expression, color = Expression)) +
          geom_point(alpha = 0.8) +
          scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) +
          theme_publication() + labs(x = "Sample", y = "Gene") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        plot.new(); text(0.5, 0.5, "No genes available", cex = 1.2)
      }
      
    } else if (input$custom_plot_type == "violin") {
      present <- intersect(rownames(logCPM), genes)
      if (length(present) == 0) {
        plot.new()
        text(0.5, 0.5, paste("None of these genes found:\n", paste(genes, collapse = ", ")), cex = 1.2)
        return()
      }
      
      grp_info <- analysisResults()$y$samples
      grp <- grp_info$group
      names(grp) <- rownames(grp_info)
      
      y_label <- NULL
      mat <- NULL
      if (identical(input$sv_metric, "per-sample logFC vs reference")) {
        ref_grp <- comparisonGroups()$ref
        if (is.na(ref_grp) || !(ref_grp %in% levels(grp))) {
          plot.new(); text(0.5, 0.5, "Select a comparison first", cex = 1.2)
          return()
        }
        ref_samples <- names(grp)[grp == ref_grp]
        if (length(ref_samples) == 0) {
          plot.new(); text(0.5, 0.5, "No samples in reference group", cex = 1.2)
          return()
        }
        ref_means <- rowMeans(logCPM[present, ref_samples, drop = FALSE])
        mat <- sweep(logCPM[present, , drop = FALSE], 1, ref_means, "-")
        y_label <- paste0("logFC vs ", ref_grp)
      } else {
        mat <- logCPM[present, , drop = FALSE]
        y_label <- "Expression (logCPM)"
      }
      
      df <- as.data.frame(mat)
      df$Gene <- rownames(df)
      df_long <- pivot_longer(df, cols = -Gene, names_to = "Sample", values_to = "Y")
      df_long$Group <- factor(grp[df_long$Sample], levels = levels(grp))
      
      ggplot(df_long, aes(x = Group, y = Y, fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, color = "black") +
        geom_jitter(width = 0.1, size = 1, alpha = 0.6, color = "black") +
        facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
        scale_fill_brewer(palette = "Set2") +
        theme_publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 13, face = "bold"),
              legend.position = "none",
              panel.spacing.y = unit(0.8, "lines")) +
        labs(x = "Group", y = y_label)
    }
  }, height = function(){
    if (input$custom_plot_type == "violin") {
      genes <- selectedGenes()
      if (length(genes) == 0) return(400)
      logCPM_check <- try(analysisResults()$logCPM, silent = TRUE)
      if (inherits(logCPM_check, "try-error")) return(400)
      n <- max(1, length(intersect(genes, rownames(logCPM_check))))
      max(400, min(6000, 280 * n))
    } else {
      600
    }
  })
  
  # Download handlers
  output$download_raw_counts <- downloadHandler(
    filename = function() { paste0("raw_counts_matrix_", Sys.Date(), ".csv") },
    content = function(file) {
      mat <- analysisResults()$counts_raw
      df <- as.data.frame(mat)
      df <- cbind(Gene = rownames(df), df)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$download_results <- downloadHandler(
    filename = "DE_results_all_comparisons.xlsx",
    content = function(file) {
      wb <- createWorkbook()
      for(comp in names(analysisResults()$comparisons)){
        addWorksheet(wb, comp)
        writeData(wb, comp, analysisResults()$comparisons[[comp]])
      }
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_siggenes <- downloadHandler(
    filename = function() { paste0("sig_genes_", gsub("/", "-", input$selected_comparison), ".csv") },
    content = function(file) {
      comp_df <- analysisResults()$comparisons[[input$selected_comparison]]
      logfc_thresh <- max(0.01, input$logfc_threshold)
      sig_df <- comp_df %>% filter(FDR <= input$fdr_threshold & abs(logFC) >= logfc_thresh)
      write.csv(sig_df, file, row.names = TRUE)
    }
  )
  
  output$download_gsea_results <- downloadHandler(
    filename = function() { 
      collection_name <- gsub(":", "_", input$msigdb_collection)
      paste0("GSEA_", collection_name, "_", Sys.Date(), ".xlsx") 
    },
    content = function(file) {
      req(analysisResults())
      fgsea_list <- analysisResults()$fgsea_list
      has_data <- any(sapply(fgsea_list, function(x) !is.null(x) && nrow(x) > 0))
      
      if (!has_data) {
        wb <- createWorkbook()
        addWorksheet(wb, "Info")
        writeData(wb, "Info", data.frame(
          Message = "No significant GSEA pathways found for any comparison.",
          Collection = input$msigdb_collection,
          Note = "Try adjusting FDR threshold or selecting a different gene set collection."
        ))
        saveWorkbook(wb, file, overwrite = TRUE)
        return()
      }
      
      wb <- createWorkbook()
      for(comp in names(fgsea_list)){
        fgsea_res <- fgsea_list[[comp]]
        if (!is.null(fgsea_res) && nrow(fgsea_res) > 0) {
          fgsea_export <- fgsea_res %>%
            mutate(
              leadingEdge = sapply(leadingEdge, function(x) {
                if (length(x) > 0) paste(x, collapse = "; ") else ""
              }),
              pval = round(pval, 6),
              padj = round(padj, 6),
              ES = round(ES, 4),
              NES = round(NES, 4)
            ) %>% arrange(padj)
          
          sheet_name <- substr(comp, 1, 31)
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, fgsea_export)
          addFilter(wb, sheet_name, row = 1, cols = 1:ncol(fgsea_export))
        }
      }
      
      summary_data <- data.frame(
        Comparison = names(fgsea_list),
        SignificantPathways = sapply(fgsea_list, function(x) {
          if (is.null(x)) return(0)
          sum(x$padj < 0.05, na.rm = TRUE)
        }),
        TotalPathwaysTested = sapply(fgsea_list, function(x) {
          if (is.null(x)) return(0)
          nrow(x)
        })
      )
      
      addWorksheet(wb, "Summary", gridLines = TRUE)
      writeData(wb, "Summary", summary_data)
      writeData(wb, "Summary", 
                data.frame(Info = c(
                  paste("Gene Set Collection:", input$msigdb_collection),
                  paste("Date:", Sys.Date()),
                  paste("FDR Threshold:", input$fdr_threshold),
                  paste("LogFC Threshold:", input$logfc_threshold),
                  paste("Batch Corrected:", analysisResults()$batch_corrected)
                )), 
                startRow = nrow(summary_data) + 3)
      
      worksheetOrder(wb) <- c(length(names(wb)), 1:(length(names(wb))-1))
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  output$download_qc <- downloadHandler(
    filename = function() { paste0("qc_", input$qc_plot_type, "_", Sys.Date(), ".png") },
    content = function(file) {
      req(analysisResults())
      
      if (input$qc_plot_type == "mito") {
        mito_df <- analysisResults()$mito_fraction_df
        if (analysisResults()$has_multiple_batches) {
          plt <- ggplot(mito_df, aes(x = Sample, y = MitoFraction, fill = Batch)) +
            scale_fill_brewer(palette = "Set2")
        } else {
          plt <- ggplot(mito_df, aes(x = Sample, y = MitoFraction)) +
            scale_fill_manual(values = APP_COLORS$primary)
        }
        plt <- plt +
          geom_bar(stat = "identity", alpha = 0.8) +
          geom_hline(yintercept = 0.1, linetype = "dashed", color = APP_COLORS$danger) +
          ylim(0, 1) +
          theme_publication() + labs(x = "Sample", y = "Mitochondrial Fraction") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(file, plt, width = 10, height = 6, dpi = 300)
        
      } else if (input$qc_plot_type == "boxplot") {
        df <- analysisResults()$norm_counts_long
        if (analysisResults()$has_multiple_batches) {
          plt <- ggplot(df, aes(x = Sample, y = Log2_CPM, fill = Batch)) +
            geom_boxplot(alpha = 0.6, outlier.color = APP_COLORS$danger) +
            scale_fill_brewer(palette = "Set2")
        } else {
          plt <- ggplot(df, aes(x = Sample, y = Log2_CPM)) +
            geom_boxplot(alpha = 0.6, outlier.color = APP_COLORS$danger, fill = APP_COLORS$primary)
        }
        plt <- plt +
          theme_publication() + labs(x = "Sample", y = "Log2 CPM") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(file, plt, width = 10, height = 6, dpi = 300)
        
      } else if (input$qc_plot_type == "violin") {
        df <- analysisResults()$norm_counts_long
        if (analysisResults()$has_multiple_batches) {
          plt <- ggplot(df, aes(x = Sample, y = Log2_CPM, fill = Batch)) +
            geom_violin(alpha = 0.6) +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA) +
            scale_fill_brewer(palette = "Set2")
        } else {
          plt <- ggplot(df, aes(x = Sample, y = Log2_CPM)) +
            geom_violin(alpha = 0.6, fill = APP_COLORS$primary) +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.shape = NA)
        }
        plt <- plt +
          theme_publication() + labs(x = "Sample", y = "Log2 CPM") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(file, plt, width = 10, height = 6, dpi = 300)
        
      } else if (input$qc_plot_type == "pca") {
        df <- analysisResults()$pca_df
        percentVar <- analysisResults()$percentVar
        has_batches <- analysisResults()$has_multiple_batches
        
        if (has_batches) {
          plt <- ggplot(df, aes(x = PC1, y = PC2, color = group, shape = batch, label = Sample)) +
            geom_point(size = 5, alpha = 0.8) + 
            ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
            scale_color_brewer(palette = "Set2") +
            scale_shape_manual(values = c(16, 17, 15, 18, 8, 3)) +
            guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
        } else {
          plt <- ggplot(df, aes(x = PC1, y = PC2, color = group, label = Sample)) +
            geom_point(size = 5, alpha = 0.8) + 
            ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
            scale_color_brewer(palette = "Set2")
        }
        
        plt <- plt +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
          theme_publication() +
          theme(legend.title = element_blank(), legend.position = "top")
        
        if (input$show_ellipses) {
          plt <- plt + stat_ellipse(aes(color = group), type = "t", level = 0.95, linewidth = 1)
        }
        ggsave(file, plt, width = 10, height = 8, dpi = 300)
      }
    }
  )
  
  output$download_de <- downloadHandler(
    filename = function() { paste0("de_", input$de_plot_type, "_", gsub("/", "-", input$selected_comparison), "_", Sys.Date(), ".png") },
    content = function(file) {
      req(analysisResults())
      
      if (input$de_plot_type == "volcano") {
        df <- selectedDE()
        pcol <- input$volcano_pval_type
        pthresh <- if (pcol == "FDR") input$fdr_threshold else input$pval_threshold
        logfc_thresh <- max(0.01, input$logfc_threshold)
        df$Significance <- ifelse(df[[pcol]] <= pthresh & abs(df$logFC) >= logfc_thresh, "Significant", "Not Significant")
        df$yval <- -log10(df[[pcol]])
        df <- df[is.finite(df$yval) & is.finite(df$logFC), ]
        sig_df <- df[df$Significance == "Significant", ]
        labels_df <- data.frame()
        if (nrow(sig_df) > 60) {
          left_labels <- head(sig_df %>% filter(logFC < 0) %>% arrange(!!rlang::sym(pcol), desc(abs(logFC))), 30)
          right_labels <- head(sig_df %>% filter(logFC > 0) %>% arrange(!!rlang::sym(pcol), desc(abs(logFC))), 30)
          labels_df <- bind_rows(left_labels, right_labels)
        } else {
          labels_df <- sig_df
        }
        plt <- ggplot(df, aes(x = logFC, y = yval, color = Significance)) +
          geom_point(alpha = 0.6, size = 2) +
          scale_color_manual(values = c("Significant" = APP_COLORS$danger, "Not Significant" = "black"))
        if (nrow(labels_df) > 0) {
          plt <- plt + ggrepel::geom_text_repel(data = labels_df, aes(label = Gene), size = 3, max.overlaps = 100, color = "black")
        }
        plt <- plt +
          geom_hline(yintercept = -log10(pthresh), linetype = "dashed", color = "black") +
          {if(logfc_thresh > 0) geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") else NULL} +
          coord_cartesian(ylim = if (is.na(input$volcano_ylim_max)) NULL else c(0, input$volcano_ylim_max)) +
          labs(title = input$selected_comparison, x = "Log2 Fold Change", y = paste0("-Log10(", pcol, ")")) +
          theme_publication_volcano_ma() + theme(legend.position = "top")
        ggsave(file, plt, width = 10, height = 8, dpi = 300)
        
      } else if (input$de_plot_type == "ma") {
        df <- selectedDE()
        logfc_thresh <- max(0.01, input$logfc_threshold)
        df$Significance <- ifelse(df$FDR <= input$fdr_threshold & abs(df$logFC) >= logfc_thresh, "Significant", "Not Significant")
        avg_logCPM <- rowMeans(analysisResults()$logCPM)
        df$avg_logCPM <- avg_logCPM[match(df$Gene, rownames(analysisResults()$logCPM))]
        plt <- ggplot(df, aes(x = avg_logCPM, y = logFC, color = Significance)) +
          geom_point(alpha = 0.6, size = 2) +
          scale_color_manual(values = c("Significant" = APP_COLORS$danger, "Not Significant" = "black")) +
          {if(logfc_thresh > 0) geom_hline(yintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "black") else NULL} +
          theme_publication_volcano_ma() + 
          labs(title = input$selected_comparison, x = "Average Log2 CPM", y = "Log2 Fold Change") +
          theme(legend.position = "top")
        ggsave(file, plt, width = 10, height = 8, dpi = 300)
        
      } else if (input$de_plot_type == "heatmap") {
        logCPM <- analysisResults()$logCPM
        mode <- input$sig_heatmap_mode
        sig_genes <- if (identical(mode, "union")) analysisResults()$sig_genes_union else analysisResults()$sig_genes_by_comp[[input$selected_comparison]]
        sig_logCPM <- logCPM[rownames(logCPM) %in% sig_genes, , drop = FALSE]
        if(nrow(sig_logCPM) > 1){
          png(file, width = 1400, height = 1200, res = 150)
          if (analysisResults()$has_multiple_batches) {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group, Batch = analysisResults()$batch)
          } else {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
          }
          rownames(annotation_col) <- colnames(sig_logCPM)
          main_title <- if (identical(mode, "union")) "Significant Genes (All)" else paste("Sig Genes:", input$selected_comparison)
          pheatmap(as.matrix(sig_logCPM), annotation_col = annotation_col, scale = "row", main = main_title,
                   fontsize_row = max(6, 10 - (nrow(sig_logCPM)/20)), fontsize_col = 11,
                   show_rownames = auto_show_rownames(nrow(sig_logCPM)),
                   color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
          dev.off()
        }
      }
    }
  )
  
  output$download_ora <- downloadHandler(
    filename = function() { paste0("ora_", ifelse(input$ora_plot_type == "dot", "dotplot", "barplot"), ".png") },
    content = function(file) {
      enrich_results <- analysisResults()$enrich_list[[input$selected_comparison]]
      plots <- list()
      show_up <- input$ora_direction %in% c("both", "up")
      show_down <- input$ora_direction %in% c("both", "down")
      
      if(show_up && !is.null(enrich_results$up) && nrow(as.data.frame(enrich_results$up)) > 0){
        if (input$ora_plot_type == "dot") {
          p_up <- dotplot(enrich_results$up, showCategory = 20) + ggtitle(paste("Upregulated:", input$selected_comparison)) +
            scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + theme_publication()
        } else {
          p_up <- barplot(enrich_results$up, showCategory = 20) + ggtitle(paste("Upregulated:", input$selected_comparison)) +
            scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + theme_publication()
        }
        plots[["up"]] <- p_up
      }
      
      if(show_down && !is.null(enrich_results$down) && nrow(as.data.frame(enrich_results$down)) > 0){
        if (input$ora_plot_type == "dot") {
          p_down <- dotplot(enrich_results$down, showCategory = 20) + ggtitle(paste("Downregulated:", input$selected_comparison)) +
            scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + theme_publication()
        } else {
          p_down <- barplot(enrich_results$down, showCategory = 20) + ggtitle(paste("Downregulated:", input$selected_comparison)) +
            scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + theme_publication()
        }
        plots[["down"]] <- p_down
      }
      
      if(length(plots) == 1){
        ggsave(file, plots[[1]], width = 12, height = 10, dpi = 300)
      } else if(length(plots) > 1){
        combined <- grid.arrange(plots[["up"]], plots[["down"]], ncol = 1)
        save_grob_png(combined, file, width = 1400, height = 1600, res = 150)
      }
    }
  )
  
  output$download_gsea <- downloadHandler(
    filename = function() { paste0("gsea_", ifelse(input$gsea_plot_type == "dot", "dotplot", "barplot"), ".png") },
    content = function(file) {
      fgseaRes <- analysisResults()$fgsea_list[[input$selected_comparison]]
      if (is.null(fgseaRes) || nrow(fgseaRes) == 0) return()
      topN <- max(1, as.integer(input$fgsea_topN))
      fg_up <- fgseaRes %>% filter(padj < 0.05, NES > 0) %>% arrange(padj) %>% head(topN)
      fg_down <- fgseaRes %>% filter(padj < 0.05, NES < 0) %>% arrange(padj) %>% head(topN)
      topPathways <- bind_rows(fg_up, fg_down)
      if (nrow(topPathways) > 0) {
        if (input$gsea_plot_type == "dot") {
          plt <- ggplot(topPathways, aes(x = NES, y = reorder(pathway, NES), size = size, color = -log10(padj))) +
            geom_point(alpha = 0.8) + scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) +
            theme_publication() + labs(title = paste("GSEA:", input$selected_comparison), x = "NES", y = "Pathway") +
            theme(axis.text.y = element_text(size = 11, face = "bold"))
        } else {
          plt <- ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
            geom_col(alpha = 0.8) + coord_flip() +
            labs(title = paste("GSEA:", input$selected_comparison), x = "Pathway", y = "NES") +
            scale_fill_gradient(low = APP_COLORS$down, high = APP_COLORS$up) + theme_publication() +
            theme(axis.text.y = element_text(size = 11, face = "bold"))
        }
        ggsave(file, plt, width = 12, height = 10, dpi = 300)
      }
    }
  )
  
  output$download_pathway_view <- downloadHandler(
    filename = function() { paste0("pathway_", input$pathway_view, ".png") },
    content = function(file) {
      req(analysisResults(), input$selected_pathways)
      
      if (input$pathway_view == "barplot") {
        mode <- input$pathway_barplot_mode
        fgsea_list <- analysisResults()$fgsea_list
        comp_names <- names(fgsea_list)
        if (mode == "selected") comp_names <- input$selected_comparison
        else if (mode == "vs_control") {
          grp_levels <- levels(analysisResults()$y$samples$group)
          control_group <- grp_levels[1]
          comp_names <- grep(paste0("^", control_group, "_vs_|_vs_", control_group, "$"), comp_names, value = TRUE)
        }
        pathway_data <- lapply(comp_names, function(comp) {
          fgsea_res <- fgsea_list[[comp]]
          if (is.null(fgsea_res)) return(NULL)
          parts <- strsplit(comp, "_vs_")[[1]]
          fgsea_res %>% filter(pathway %in% input$selected_pathways) %>%
            mutate(Comparison = comp, Reference = parts[1], Group = parts[2],
                   sig_label = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))
        })
        pathway_data <- bind_rows(pathway_data)
        if (nrow(pathway_data) > 0) {
          p <- ggplot(pathway_data, aes(x = pathway, y = NES, fill = Group)) +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.8) + 
            geom_col(position = position_dodge(width = 0.8), color = NA, width = 0.75, alpha = 0.8) +
            geom_text(aes(label = sig_label, y = NES + sign(NES) * max(abs(NES)) * 0.08),
                      position = position_dodge(width = 0.8), size = 5, fontface = "bold") +
            coord_flip() + scale_fill_brewer(palette = "Set2") + theme_publication() + 
            labs(x = "Pathway", y = "Normalized Enrichment Score (NES)")
          if (mode == "all_pairwise") p <- p + facet_wrap(~ paste0("Ref: ", Reference), ncol = 2)
          ggsave(file, p, width = 12, height = 10, dpi = 300)
        }
      } else if (input$pathway_view == "traces") {
        comp_df <- analysisResults()$comparisons[[input$selected_comparison]]
        pvals <- comp_df$PValue
        min_nonzero <- min(pvals[pvals > 0 & is.finite(pvals)], na.rm = TRUE)
        pvals[!is.finite(pvals) | pvals <= 0] <- min_nonzero
        stats <- sign(comp_df$logFC) * -log10(pvals) * abs(comp_df$logFC)
        names(stats) <- comp_df$Gene
        stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
        plots <- lapply(input$selected_pathways, function(pw) {
          if (pw %in% names(analysisResults()$pathways_list)) {
            plotEnrichment(analysisResults()$pathways_list[[pw]], stats) + labs(title = pw) + theme_publication()
          } else NULL
        })
        plots <- Filter(Negate(is.null), plots)
        if(length(plots) == 1) ggsave(file, plots[[1]], width = 10, height = 8, dpi = 300)
        else if(length(plots) > 1){
          combined <- do.call(grid.arrange, c(plots, ncol = 1))
          save_grob_png(combined, file, width = 1400, height = 400 * length(plots), res = 150)
        }
      } else if (input$pathway_view == "heatmap") {
        logCPM <- analysisResults()$logCPM
        pathways_list <- analysisResults()$pathways_list
        gene_to_pathway <- list()
        for (pw in input$selected_pathways) {
          if (pw %in% names(pathways_list)) {
            genes <- pathways_list[[pw]]
            for (g in genes) {
              if (g %in% rownames(logCPM)) {
                if (is.null(gene_to_pathway[[g]])) gene_to_pathway[[g]] <- pw
                else gene_to_pathway[[g]] <- paste(gene_to_pathway[[g]], pw, sep = "; ")
              }
            }
          }
        }
        selected_genes <- names(gene_to_pathway)
        selected_logCPM <- logCPM[selected_genes, , drop = FALSE]
        if(nrow(selected_logCPM) > 1){
          png(file, width = 1400, height = 1200, res = 150)
          if (analysisResults()$has_multiple_batches) {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group, Batch = analysisResults()$batch)
          } else {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
          }
          rownames(annotation_col) <- colnames(selected_logCPM)
          annotation_row <- data.frame(Pathway = unlist(gene_to_pathway[rownames(selected_logCPM)]))
          rownames(annotation_row) <- rownames(selected_logCPM)
          pheatmap(as.matrix(selected_logCPM), annotation_col = annotation_col, annotation_row = annotation_row,
                   scale = "row", show_rownames = auto_show_rownames(nrow(selected_logCPM)),
                   color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
          dev.off()
        }
      }
    }
  )
  
  output$download_custom <- downloadHandler(
    filename = function() { paste0("custom_genes_", input$custom_plot_type, ".png") },
    content = function(file) {
      req(analysisResults())
      genes <- selectedGenes()
      if (length(genes) == 0) return()
      logCPM <- analysisResults()$logCPM
      
      if (input$custom_plot_type == "heatmap") {
        selected_logCPM <- logCPM[intersect(rownames(logCPM), genes), , drop = FALSE]
        if(nrow(selected_logCPM) > 1){
          png(file, width = 1400, height = 1200, res = 150)
          if (analysisResults()$has_multiple_batches) {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group, Batch = analysisResults()$batch)
          } else {
            annotation_col <- data.frame(Group = analysisResults()$y$samples$group)
          }
          rownames(annotation_col) <- colnames(selected_logCPM)
          pheatmap(as.matrix(selected_logCPM), annotation_col = annotation_col, scale = "row",
                   show_rownames = auto_show_rownames(nrow(selected_logCPM)),
                   color = colorRampPalette(c(APP_COLORS$down, "white", APP_COLORS$up))(100))
          dev.off()
        }
      } else if (input$custom_plot_type == "dotplot") {
        dot_logCPM <- logCPM[intersect(rownames(logCPM), genes), , drop = FALSE]
        if(nrow(dot_logCPM) > 0){
          dot_df <- as.data.frame(dot_logCPM)
          dot_df$Gene <- rownames(dot_df)
          dot_df_long <- pivot_longer(dot_df, cols = -Gene, names_to = "Sample", values_to = "Expression")
          plt <- ggplot(dot_df_long, aes(x = Sample, y = Gene, size = Expression, color = Expression)) +
            geom_point(alpha = 0.8) + scale_color_gradient(low = APP_COLORS$down, high = APP_COLORS$up) +
            theme_publication() + labs(x = "Sample", y = "Gene") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          ggsave(file, plt, width = 10, height = 8, dpi = 300)
        }
      } else if (input$custom_plot_type == "violin") {
        present <- intersect(rownames(logCPM), genes)
        if (length(present) == 0) return()
        grp_info <- analysisResults()$y$samples
        grp <- grp_info$group
        names(grp) <- rownames(grp_info)
        y_label <- NULL
        mat <- NULL
        if (identical(input$sv_metric, "per-sample logFC vs reference")) {
          ref_grp <- comparisonGroups()$ref
          ref_samples <- names(grp)[grp == ref_grp]
          ref_means <- rowMeans(logCPM[present, ref_samples, drop = FALSE])
          mat <- sweep(logCPM[present, , drop = FALSE], 1, ref_means, "-")
          y_label <- paste0("logFC vs ", ref_grp)
        } else {
          mat <- logCPM[present, , drop = FALSE]
          y_label <- "Expression (logCPM)"
        }
        df <- as.data.frame(mat)
        df$Gene <- rownames(df)
        df_long <- pivot_longer(df, cols = -Gene, names_to = "Sample", values_to = "Y")
        df_long$Group <- factor(grp[df_long$Sample], levels = levels(grp))
        plt <- ggplot(df_long, aes(x = Group, y = Y, fill = Group)) +
          geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
          geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8, color = "black") +
          geom_jitter(width = 0.1, size = 1, alpha = 0.6, color = "black") +
          facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
          scale_fill_brewer(palette = "Set2") + theme_publication() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
          labs(x = "Group", y = y_label)
        h_px <- max(500, min(6000, 280 * length(present)))
        ggsave(file, plt, width = 10, height = h_px / 100, dpi = 300)
      }
    }
  )
}

# Run the application 
shinyApp(ui, server)