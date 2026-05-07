# ============================================================
# Shiny RNA-seq / Omics Explorer (Enhanced version)
# Developed by Janan Gawra
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(DT)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  library(ggrepel)
  library(readxl)
  library(colourpicker)
  library(ComplexHeatmap)
  library(circlize)
})
options(shiny.maxRequestSize = 2000 * 1024^2)

# ---- Helper functions ----
read_table_safely <- function(path, header = TRUE) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  df <- tryCatch({
    if (ext %in% c("xlsx", "xls")) read_excel(path)
    else if (ext %in% c("tsv", "txt")) read.delim(path, header = header, check.names = FALSE)
    else if (ext %in% c("csv")) read.csv(path, header = header, check.names = FALSE)
    else stop("Unsupported file type")
  }, error = function(e) NULL)
  if (inherits(df, "tbl_df")) df <- as.data.frame(df)
  df[] <- lapply(df, function(x) if (is.numeric(x)) x else as.character(x))
  df
}

scale_rows <- function(mat) t(scale(t(mat)))

# ---- UI ----
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # Title + footer
  fluidRow(
    column(
      12,
      div(
        h1("Shiny RNA-seq / Omics Explorer", style = "font-weight:bold;"),
        div(
          HTML(
            paste(
              "<span style='color:blue; font-size:0.9em;'>Developed by Janan Gawra</span>",
              "<a href='https://www.linkedin.com/in/janangawra/' target='_blank' style='font-size:0.9em;'> | <i class='fab fa-linkedin'></i> LinkedIn Profile</a>"
            )
          ),
          style = "text-align:left; margin-bottom:10px;"
        )
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      h4("1) Upload files"),
      fileInput("counts_file", "Counts / normalized expression",
                accept = c('.csv', '.tsv', '.txt', '.xlsx', '.xls')),
      checkboxInput("counts_header", "First row has headers", TRUE),
      uiOutput("counts_gene_col_ui"),
      hr(),
      
      fileInput("coldata_file", "Sample metadata (colData)",
                accept = c('.csv', '.tsv', '.txt', '.xlsx', '.xls')),
      checkboxInput("coldata_header", "First row has headers", TRUE),
      uiOutput("coldata_sample_id_ui"),
      hr(),
      
      h4("2) Coloring options"),
      uiOutput("annot_factor_ui"),
      uiOutput("color_picker_ui"),
      selectInput("highlight_sample", "Highlight sample", choices = "None"),
      hr(),
      
      h4("3) Plot customization"),
      textInput("heatmap_title", "Heatmap title", value = "Heatmap of Gene Expression"),
      textInput("pca_title", "PCA title", value = "Principal Component Analysis"),
      textInput("pca_xlab", "PCA X-axis label", value = "PC1"),
      textInput("pca_ylab", "PCA Y-axis label", value = "PC2"),
      textInput("geneplot_title", "Gene Expression title", value = "Expression per Group"),
      textInput("geneplot_ylab", "Gene Expression Y-axis label", value = "Expression"),
      hr(),
      checkboxInput("log_transform", "Log2(x + 1) transform", TRUE)
    ),
    
    mainPanel(
      tabsetPanel(
        # ---- DESCRIPTION TAB ----
        tabPanel("📘 Description / Help",
                 div(
                   h3("📘 How to Use the Shiny RNA-seq / Omics Explorer"),
                   p("This app lets you explore RNA-seq or other omics data interactively, creating publication-ready heatmaps, PCA plots, and gene expression boxplots."),
                   hr(),
                   h4("1️⃣ Uploading Files"),
                   tags$ul(
                     tags$li(strong("Counts / normalized expression file:"), 
                             "Upload a CSV, TSV, TXT, or Excel file (genes × samples)."),
                     tags$li(strong("Gene ID column:"), "Select the column containing gene identifiers."),
                     tags$li(strong("Sample metadata (colData):"), 
                             "Upload metadata with one row per sample (conditions, replicates, etc.)."),
                     tags$li(strong("Sample ID column:"), 
                             "Select the column that uniquely identifies each sample (must match column names in the count matrix).")
                   ),
                   hr(),
                   h4("2️⃣ Coloring Options"),
                   tags$ul(
                     tags$li("Select a metadata factor to color plots."),
                     tags$li("Manually pick colors for each group."),
                     tags$li("Optionally highlight a single sample.")
                   ),
                   hr(),
                   h4("3️⃣ Heatmap Tab"),
                   tags$ul(
                     tags$li("Choose all genes, top variable genes, or a custom list."),
                     tags$li("Select scaling: none, Z-score rows, or columns."),
                     tags$li("Enable/disable clustering for rows/columns."),
                     tags$li("Click 'Render heatmap' to visualize.")
                   ),
                   hr(),
                   h4("4️⃣ PCA Tab"),
                   tags$ul(
                     tags$li("Performs PCA on log2-transformed and optionally scaled data."),
                     tags$li("Colors samples by the selected metadata factor."),
                     tags$li("Edit axis labels and titles.")
                   ),
                   hr(),
                   h4("5️⃣ Gene Expression Tab"),
                   tags$ul(
                     tags$li("Type or select a gene ID to view expression."),
                     tags$li("Choose plot type: Boxplot, Violin, or Jitter."),
                     tags$li("Log-transform and show points optionally.")
                   ),
                   hr(),
                   h4("📂 Example Input Format"),
                   pre("
Counts file:
Gene_ID, CD_CTRL_1, CD_CTRL_2, CD_C1_1, CD_C1_2
ENSDARG0000000018, 123, 245, 312, 221
ENSDARG0000000321, 432, 534, 312, 298
...

colData file:
SampleName, Condition, Replicate
CD_CTRL_1, Control, R1
CD_C1_1, Low_dose, R1
CD_C3_1, Medium_dose, R1
CD_C5_1, High_dose, R1
"),
                   hr(),
                   h4("❓ Need Help?"),
                   HTML("<p>For questions or feedback, contact <b>Janan Gawra</b> on 
                        <a href='https://www.linkedin.com/in/janangawra/' target='_blank'>LinkedIn</a>.</p>")
                 )
        ),
        
        # ---- EXISTING TABS ----
        tabPanel("Data preview",
                 fluidRow(
                   column(6, h4("Counts preview"), DTOutput("counts_tbl")),
                   column(6, h4("colData preview"), DTOutput("col_tbl"))
                 ),
                 hr(),
                 verbatimTextOutput("match_summary")
        ),
        
        tabPanel("Heatmap",
                 fluidRow(
                   column(3,
                          selectInput("heatmap_gene_source", "Gene subset",
                                      c("All genes"="all","Top variable genes"="hvg","Custom list"="custom")),
                          conditionalPanel("input.heatmap_gene_source == 'hvg'",
                                           numericInput("hvg_n", "Top variable genes", 500, min=10)),
                          conditionalPanel("input.heatmap_gene_source == 'custom'",
                                           textAreaInput("custom_genes", "Genes (comma-separated)", rows=3)),
                          selectInput(
                            "scale_option",
                            "Scaling option",
                            choices = c("None" = "none", "Z-score rows" = "rows", "Z-score columns" = "cols"),
                            selected = "rows"
                          ),
                          checkboxInput("cluster_rows", "Cluster rows", TRUE),
                          checkboxInput("cluster_cols", "Cluster columns", TRUE),
                          selectInput("heatmap_palette", "Palette",
                                      c("RdBu","PuOr","PiYG","Spectral","RdYlBu")),
                          checkboxInput("show_rownames", "Show gene labels", FALSE),
                          checkboxInput("show_colnames", "Show sample labels", TRUE),
                          actionButton("plot_heatmap", "Render heatmap", class = "btn-primary")
                   ),
                   column(9, plotOutput("heatmap_plot", height="700px"))
                 )
        ),
        
        tabPanel("PCA",
                 fluidRow(
                   column(3,
                          checkboxInput("pca_scale", "Scale genes", TRUE),
                          selectInput("pca_label", "Point labels", c("None","SampleID")),
                          numericInput("pca_pc_x", "X-axis PC", 1, min=1),
                          numericInput("pca_pc_y", "Y-axis PC", 2, min=1)
                   ),
                   column(9, plotlyOutput("pca_plot", height="650px"),
                          verbatimTextOutput("pca_varexp"))
                 )
        ),
        
        tabPanel("Gene expression",
                 fluidRow(
                   column(3,
                          uiOutput("gene_select_ui"),
                          selectInput("expr_geom", "Plot type", c("Boxplot","Violin","Jitter")),
                          checkboxInput("expr_log", "Log2(x + 1) transform", TRUE),
                          checkboxInput("expr_points", "Show points", TRUE)
                   ),
                   column(9, plotlyOutput("gene_plot", height="650px"))
                 )
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  # Load counts
  counts_raw <- reactive({
    inFile <- input$counts_file; req(inFile)
    read_table_safely(inFile$datapath, header = input$counts_header)
  })
  output$counts_gene_col_ui <- renderUI({
    df <- counts_raw(); req(df)
    selectInput("counts_gene_col", "Gene ID column", choices = colnames(df), selected = colnames(df)[1])
  })
  
  counts <- reactive({
    df <- counts_raw(); req(input$counts_gene_col)
    rn <- make.unique(as.character(df[[input$counts_gene_col]]))
    df[[input$counts_gene_col]] <- NULL
    df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
    df <- df[, vapply(df, function(x) all(is.finite(x)|is.na(x)), logical(1)), drop=FALSE]
    rownames(df) <- rn
    df
  })
  
  # Load colData
  coldata_raw <- reactive({
    inFile <- input$coldata_file; req(inFile)
    read_table_safely(inFile$datapath, header = input$coldata_header)
  })
  output$coldata_sample_id_ui <- renderUI({
    df <- coldata_raw()
    uniq_cols <- names(which(vapply(df, function(x) length(unique(x)) == nrow(df), logical(1))))
    default_col <- if (length(uniq_cols) > 0) uniq_cols[1] else colnames(df)[1]
    selectInput("col_sample_col", "Sample ID column",
                choices = colnames(df), selected = default_col)
  })
  
  matched <- reactive({
    cts <- counts(); cd <- coldata_raw(); req(input$col_sample_col)
    sample_ids <- as.character(cd[[input$col_sample_col]])
    validate(need(anyDuplicated(sample_ids) == 0, "Sample IDs must be unique."))
    common <- intersect(colnames(cts), sample_ids)
    validate(need(length(common) > 0, "No overlapping sample IDs found."))
    cd_match <- cd[match(common, sample_ids), , drop = FALSE]
    rownames(cd_match) <- common
    cts_match <- cts[, common, drop = FALSE]
    for (col in colnames(cd_match)) {
      if (is.character(cd_match[[col]]) || is.factor(cd_match[[col]]))
        cd_match[[col]] <- factor(cd_match[[col]], levels = unique(cd[[col]]))
    }
    list(counts = cts_match, coldata = cd_match)
  })
  # ---- Data Preview ----
  output$counts_tbl <- renderDT({
    req(counts_raw())
    datatable(head(counts_raw(), 50), options = list(scrollX = TRUE))
  })
  
  output$col_tbl <- renderDT({
    req(coldata_raw())
    datatable(head(coldata_raw(), 50), options = list(scrollX = TRUE))
  })
  
  output$match_summary <- renderPrint({
    m <- matched()
    cat("✅ Data successfully loaded!\n")
    cat("Genes:", nrow(m$counts), "\nSamples:", ncol(m$counts), "\n")
  })
  
  # Preview
  output$counts_tbl <- renderDT({datatable(head(counts_raw(),50), options=list(scrollX=TRUE))})
  output$col_tbl <- renderDT({datatable(head(coldata_raw(),50), options=list(scrollX=TRUE))})
  output$match_summary <- renderPrint({
    m <- matched(); cat("Genes:",nrow(m$counts),"\nSamples:",ncol(m$counts),"\n")
  })
  
  # Color options
  output$annot_factor_ui <- renderUI({
    m <- matched()
    selectInput("annot_factor", "Color by factor", choices = c("None", colnames(m$coldata)), selected="None")
  })
  output$color_picker_ui <- renderUI({
    req(input$annot_factor)
    if (input$annot_factor=="None") return(NULL)
    f <- as.factor(matched()$coldata[[input$annot_factor]])
    levs <- unique(f)
    lapply(levs, function(lv)
      colourInput(inputId=paste0("col_",lv), label=paste("Color for",lv), value="#CCCCCC"))
  })
  user_colors <- reactive({
    if (is.null(input$annot_factor) || input$annot_factor=="None") return(NULL)
    f <- as.factor(matched()$coldata[[input$annot_factor]])
    levs <- unique(f)
    cols <- sapply(levs, function(lv) input[[paste0("col_",lv)]])
    names(cols) <- levs; cols
  })
  observe({
    updateSelectInput(session, "highlight_sample", choices = c("None", colnames(matched()$counts)))
  })
  
  # ---- ComplexHeatmap ----
  observeEvent(input$plot_heatmap, {
    output$heatmap_plot <- renderPlot({
      m <- matched(); x <- as.matrix(m$counts)
      if (input$log_transform) x <- log2(x + 1)
      if (input$heatmap_gene_source=="hvg") {
        vars <- apply(x,1,var)
        keep <- names(sort(vars,decreasing=TRUE))[seq_len(min(length(vars),input$hvg_n))]
        x <- x[keep,,drop=FALSE]
      } else if (input$heatmap_gene_source=="custom") {
        genes <- trimws(unlist(strsplit(input$custom_genes,"[,;\\n\\t]+"))); genes <- genes[genes!=""]
        x <- x[intersect(rownames(x),genes),,drop=FALSE]
      }
      
      if (input$scale_option == "rows") x <- scale_rows(x)
      else if (input$scale_option == "cols") x <- t(scale(t(x)))
      
      pal_name <- ifelse(is.null(input$heatmap_palette) || input$heatmap_palette == "",
                         "RdBu", input$heatmap_palette)
      pal <- tryCatch(rev(brewer.pal(11, pal_name)),
                      error = function(e) rev(brewer.pal(11, "RdBu")))
      
      ha <- NULL
      if (!is.null(input$annot_factor) && input$annot_factor != "None") {
        annot_name <- input$annot_factor
        ha <- HeatmapAnnotation(
          df = m$coldata[, annot_name, drop = FALSE],
          col = setNames(list(user_colors()), annot_name)
        )
      }
      
      Heatmap(
        x,
        name = "expression",
        top_annotation = ha,
        show_row_names = input$show_rownames,
        show_column_names = input$show_colnames,
        cluster_rows = input$cluster_rows,
        cluster_columns = input$cluster_cols,
        column_title = input$heatmap_title,
        col = colorRamp2(c(-2, 0, 2), pal[c(1, 6, 11)])
      )
    })
  })
  
  # PCA
  output$pca_plot <- renderPlotly({
    m <- matched(); x <- as.matrix(m$counts)
    if (input$log_transform) x <- log2(x + 1)
    if (input$pca_scale) x <- scale_rows(x)
    pr <- prcomp(t(x), center=TRUE)
    meta <- m$coldata
    df <- as.data.frame(pr$x)
    df$SampleID <- rownames(df)
    gp <- if (input$annot_factor == "None") "All" else input$annot_factor
    df$Group <- factor(meta[rownames(df), gp], levels = unique(meta[[gp]]))
    g <- ggplot(df, aes_string(paste0("PC",input$pca_pc_x), paste0("PC",input$pca_pc_y), color="Group")) +
      geom_point(size=3) + theme_minimal(base_size=14) +
      labs(color=gp, x=input$pca_xlab, y=input$pca_ylab, title=input$pca_title)
    if (input$pca_label=="SampleID") g <- g + geom_text_repel(aes(label=SampleID))
    if (!is.null(user_colors())) g <- g + scale_color_manual(values=user_colors())
    ggplotly(g)
  })
  
  output$pca_varexp <- renderPrint({
    m <- matched(); x <- as.matrix(m$counts)
    if (input$log_transform) x <- log2(x + 1)
    if (input$pca_scale) x <- scale_rows(x)
    pr <- prcomp(t(x), center=TRUE)
    v <- pr$sdev^2; cat("Variance explained (first 10 PCs):\n")
    print(round(v[1:min(10, length(v))] / sum(v) * 100, 2))
  })
  
  # Gene expression
  output$gene_select_ui <- renderUI({
    m <- matched()
    unique_genes <- unique(gsub("\\.[0-9]+$","", rownames(m$counts)))
    selectizeInput(
      "gene_select",
      "Gene",
      choices = unique_genes,
      selected = unique_genes[1],
      options = list(
        maxOptions = 2000,
        create = TRUE,           # ✅ allows manual typing of new genes
        placeholder = 'Type or select a gene name...'
      )
    )
  })
  
  
  gene_long <- reactive({
    m <- matched(); g <- trimws(input$gene_select)
    rn <- trimws(rownames(m$counts))
    match_rows <- which(g==rn | grepl(paste0("^",g,"(\\.|$)"), rn))
    validate(need(length(match_rows)>0, paste("Gene",g,"not found")))
    x <- m$counts[match_rows,,drop=FALSE]
    if (nrow(x)>1) x <- colMeans(x,na.rm=TRUE) else x <- as.numeric(x)
    if (input$expr_log) x <- log2(x + 1)
    df <- data.frame(SampleID=colnames(m$counts), expr=as.numeric(x),
                     m$coldata, check.names=FALSE)
    df
  })
  
  output$gene_plot <- renderPlotly({
    df <- gene_long(); gp <- input$annot_factor
    if (is.null(gp) || gp=="None") gp <- colnames(matched()$coldata)[1]
    df[[gp]] <- factor(df[[gp]], levels = levels(matched()$coldata[[gp]]))
    g <- ggplot(df, aes_string(gp, "expr", color=gp)) +
      theme_minimal(base_size=14) +
      labs(x=gp, y=input$geneplot_ylab, title=paste(input$geneplot_title, "-", input$gene_select))
    if (input$expr_geom=="Boxplot") g <- g + geom_boxplot(outlier.shape=NA)
    if (input$expr_geom=="Violin") g <- g + geom_violin(trim=FALSE)
    if (input$expr_points || input$expr_geom=="Jitter") g <- g + geom_jitter(width=0.15, size=2)
    if (!is.null(user_colors())) g <- g + scale_color_manual(values=user_colors())
    ggplotly(g)
  })
}

# ---- Run App ----
shinyApp(ui, server)
