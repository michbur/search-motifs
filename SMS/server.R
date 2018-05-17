library(shiny)
library(genbankr)
library(Biostrings)
library(DT)
library(mlr)
library(dplyr)
library(PCRedux)
library(ggplot2)
library(plotly)
library(reshape2)

options(shiny.maxRequestSize=30*1024^2)

options(DT.options = list(search = list(regex = TRUE, search = 'CXC'),
                          dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 15
))

source("functions.R", local = TRUE)

shinyServer(function(input, output) {
  
  raw_input <- reactive({
    
    input_dat <- NULL
    
    if (!is.null(input[["seq_file"]]))
      input_dat <- readGenBank2(input[["seq_file"]][["datapath"]])
    # input[["use_example"]]
    # isolate({
    #   if (!is.null(input[["use_example"]]))
    #     if(input[["use_example"]] > 0)
    #       input_curves <- read.csv("example_dat.csv")
    # })
    # 
    # if(exists("input_curves")) {
    #   # here we should put some tests of the input, like to maximum amount of curves
    #   # the first column should be named cycle
    #   colnames(input_curves)[1] <- "cycle"
    #   input_curves
    # } else {
    #   NULL
    # }
    
    input_dat
  })
  
  processed_input <- reactive({
    seqs <- cds(raw_input())

    all_cds <- as.vector(mcols(seqs)[["translation"]])
    
    seqs_df <- mcols(seqs)
    data.frame(select(data.frame(ranges(seqs)), -width), 
               strand = strand(seqs),
               gene_id = seqs_df[["gene_id"]],
               gene = seqs_df[["gene"]],
               #locus_tag = seqs_df[["locus_tag"]],
               #transcript_id = seqs_df[["transcript_id"]],
               fun = as.list(seqs_df[["function"]]) %>% 
                 lapply(paste0, collapse = "; ") %>% 
                 unlist,
               protein_id = seqs_df[["protein_id"]],
               product = seqs_df[["product"]],
               note = seqs_df[["note"]],
               seq = as.vector(seqs_df[["translation"]]))
  })
  
  output[["gkb_dt"]] <- DT::renderDataTable({
    my_DT(processed_input())
  })
  
  
  output$dynamic_tabset <- renderUI({
    if(is.null(raw_input())) {
      tabsetPanel(
        tabPanel(title = "Data input",
                 fluidRow(
                   column(width = 5, 
                          fileInput('seq_file', "Submit data (.fasta or .gbk file):")
                   )
                 )
        )
      )
    } else {
      tabsetPanel(
        tabPanel(title = "Results",
                 #DT::dataTableOutput("pred_table"),
                 DT::dataTableOutput("gkb_dt"),
                 #plotOutput("prediction_plot"),
                 tags$p(HTML("<h3><A HREF=\"javascript:history.go(0)\">Start a new query</A></h3>"))
        )
        # tabPanel(title = "PCRedux data",
        #          DT::dataTableOutput("processed_input_dt")),
        # tabPanel(title = "PCRedux prediction",
        #          DT::dataTableOutput("prediction_dt")),
        # tabPanel(title = "PCRedux prediction plot",
        #          plotOutput("prediction_plot"))
      )
    }
  })
})

