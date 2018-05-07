library(shiny)
library(stringi)
library(ggplot2)

source("get_motif.R", local = TRUE)

shinyServer(function(input, output, session) {
  
  
  final_motif <- reactive(
    paste0(c("(?=", 
             sapply(1L:input[["number_of_motifs"]], function(i)
               process_single_motif(input[[paste0("mot", i)]], 
                                    input[[paste0("len_min", i)]], 
                                    input[[paste0("len_max", i)]])),
             ")"),
           collapse = "")
  )
  
  motif_pos <-reactive({
    #stri_locate_all(input[["seq"]], regex = final_motif())
    starts <- gregexpr(final_motif(), input[["seq"]], perl = TRUE)[[1]]
    ends <- 1 + rev(nchar(input[["seq"]]) - as.vector(gregexpr(final_motif(), stringi::stri_reverse(input[["seq"]]), perl = TRUE)[[1]]))
    
    lapply(1L:length(starts), function(i)
      starts[i]:ends[i])
  })
  
  v_pos <- reactive({
    unlist(motif_pos())
  })
  
  v_seq <- reactive({
    strsplit(input[["seq"]], "")[[1]]
  })
  
  output[["final_motif"]] <- renderPrint({
    final_motif()
  })
  
  
  output[["extracted_motifs"]] <- renderText({
    in_motif <- 1L:length(v_seq()) %in% v_pos()
    nice_seq <- v_seq()
    nice_seq[in_motif] <- paste0("<font color='#FF0000'><b>", nice_seq[in_motif], "</b></font>")
    paste0("<font face='Courier New,Courier,monospace'>", nice_seq, "</font>")
  })
  
  output[["human_readable_motif"]] <- renderPrint({
    all_variants <- lapply(1L:input[["number_of_motifs"]], function(i) 
      sapply(input[[paste0("len_min", i)]]:input[[paste0("len_max", i)]], function(j)
        paste0(rep(input[[paste0("mot", i)]], j), collapse = ""))
    )
    
    expand.grid(all_variants) %>% 
      apply(1, paste0, collapse = "") %>% 
      cat(sep = "\n")
  })
  
  output[["motif_position"]] <- renderPrint({
    v_pos()
  })
  
  output[["motif_position_plot"]] <- renderPlot({
    data.frame(pos = 1L:nchar(input[["seq"]]),
               aa =  strsplit(input[["seq"]], "")[[1]]) %>%
      mutate(mot = pos %in% v_pos(),
             pos_disc = floor(pos/80)) %>% 
      ggplot(aes(x = pos, y = mot, label = aa)) + 
      geom_text() +
      scale_y_discrete("Amino acids matched to the motif") +
      scale_x_continuous("Position") +
      facet_wrap(~ pos_disc, scales = "free_x", ncol = 1) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
  })
  
  # output[["color_seq"]] <- renderPrint({
  #   v_seq <- sttrsplit(input[["seq"]], "")[[""]]
  #   pos <- stri_locate_last_regex(pattern = final_motif(), input[["seq"]])[[1]]
  #   
  #   paste0("hello input is","<font color=\"#FF0000\"><b>", input$n, "</b></font>") 
  #   
  #   apply(pos, 1, function(ith_row)
  #     v_seq[ith_row[1]:ith_row[2]]
  #   )
  # })
  
  output[["motif_ui"]] <- renderUI({
    lapply(1L:input[["number_of_motifs"]], function(i) {
      list(h3(paste0("Part ", i)),
           fluidRow(
             column(4,
                    textInput(paste0("mot", i), paste0("Regular expression:"), "X")
             ),
             column(3,
                    numericInput(paste0("len_min", i), "Min. length:", value = 1, min = 1, max = 10, step = 1)
             )
             ,
             column(3,
                    numericInput(paste0("len_max", i), "Max. length:", value = 1, min = 1, max = 10, step = 1)
             )
           )
      )
    })
  })
  
  observe({
    for(i in 1L:input[["number_of_motifs"]]) {
      updateTextInput(session, paste0("mot", i), paste0("Regular expression:"), 
                      value = lapply(reactiveValuesToList(input), unclass)[[paste0("mot", i)]])
      updateNumericInput(session, paste0("len_min", i), "Min. length:", 
                         value = lapply(reactiveValuesToList(input), unclass)[[paste0("len_min", i)]])
      updateNumericInput(session, paste0("len_max", i), "Max. length:", 
                         value = lapply(reactiveValuesToList(input), unclass)[[paste0("len_max", i)]])
    }
  })
  
})
