library(shiny)

shinyUI(fluidPage(#tags$head(includeScript("ga.js")),
  #tags$style(includeCSS("./www/report.css")),
  theme = shinythemes::shinytheme("spacelab"),
  tags$style(HTML("                  
                  .shiny-input-container:not(.shiny-input-container-inline) {
                  width: 100%;
                  }
                  
                  pre{
                  background: white;
                  }
                  
                  .shiny-notification {
                  height: 100px;
                  width: 800px;
                  position:fixed;
                  top: calc(50% - 50px);;
                  left: calc(50% - 400px);;
                  }
                  ")),
  title = "Sujeet Motif Searcher",
  
  headerPanel(""),
  
  sidebarLayout(
    sidebarPanel(style = "background-color: #e0e0e0;",
                 includeMarkdown("readme.md")
  ),
  
  mainPanel(
    uiOutput("dynamic_tabset")
  )
  )
))
