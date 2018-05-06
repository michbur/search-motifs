library(shiny)
library(stringi)

source("get_motif.R", local = TRUE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  final_motif <- reactive(
    paste0(sapply(1L:input[["number_of_motifs"]], function(i)
      process_single_motif(input[[paste0("mot", i)]], 
                           input[[paste0("len", i)]][1], 
                           input[[paste0("len", i)]][2])),
      collapse = "")
  )
  
  output[["final_motif"]] <- renderPrint({
    final_motif()
  })
  
  output[["motif_position"]] <- renderPrint({
    stri_locate_all(input[["seq"]], regex = final_motif())
  })
  
  output[["color_seq"]] <- renderPrint({
    v_seq <- sttrsplit(input[["seq"]], "")[[""]]
    pos <- stri_locate_last_regex(pattern = final_motif(), input[["seq"]])[[1]]
    
    paste0("hello input is","<font color=\"#FF0000\"><b>", input$n, "</b></font>") 
    
    apply(pos, 1, function(ith_row)
      v_seq[ith_row[1]:ith_row[2]]
    )
    
  })
  
  output[["motif_ui"]] <- renderUI({
    lapply(1L:input[["number_of_motifs"]], function(i) {
      list(h3(paste0("Part ", i)),
           textInput(paste0("mot", i), paste0("Regular expression:"), "X"),
           sliderInput(paste0("len", i),
                       paste0("Length:"),
                       min = 1,
                       max = 10,
                       value = c(1, 2)))
    })
  })
  
})
