#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Motif test"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("seq", "Sequence to look for a motif", "TTNAAAANAAAANTT"),
      numericInput("number_of_motifs", "Number of motifs", value = 3, min = 1, max = 10, step = 1),
      uiOutput("motif_ui")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h3("Motif as a regular expression"),
      verbatimTextOutput("final_motif"),
      h3("Position of the match"),
      verbatimTextOutput("motif_position")
    )
  )
))
