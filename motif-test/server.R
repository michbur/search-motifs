library(shiny)
library(stringi)
library(ggplot2)

source("get_motif.R", local = TRUE)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  final_motif <- reactive(
    paste0(c("(?=", 
             sapply(1L:input[["number_of_motifs"]], function(i)
               process_single_motif(input[[paste0("mot", i)]], 
                                    input[[paste0("len", i)]][1], 
                                    input[[paste0("len", i)]][2])),
             ")"),
           collapse = "")
  )
  
  motif_pos <-reactive({
    #stri_locate_all(input[["seq"]], regex = final_motif())
    starts <- gregexpr("(?=B[A]{2,3}B)", seq, perl = TRUE)[[1]]
    ends <- rev(nchar(seq) - as.vector(gregexpr("(?=B[A]{2,3}B)", stringi::stri_reverse(seq), perl = TRUE)[[1]]))
    
    unlist(lapply(1L:length(starts), function(i)
      starts[i]:ends[i]))
  })
  
  output[["final_motif"]] <- renderPrint({
    final_motif()
  })
  
  output[["motif_position"]] <- renderPrint({
    motif_pos()
  })
  
  output[["motif_position_plot"]] <- renderPlot({
    data.frame(pos = 1L:nchar(input[["seq"]]),
               aa =  strsplit(input[["seq"]], "")[[1]]) %>%
      mutate(mot = pos %in% motif_pos()) %>% 
      ggplot(aes(x = pos, y = mot, label = aa)) + 
      geom_text()
    
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
