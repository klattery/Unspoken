library(shiny)
dir_run <- "Hello"
myout <- "con_test.csv"
work <- "home"


ui <- fluidPage(
  # theme = shinytheme("flatly"),
  tags$style("p { color: blue;}"),
  titlePanel("Unspoken Estimation Has Finished"),
  sidebarLayout(
    sidebarPanel(
      h4(p("Please click to download")),
      h3(downloadButton("downloadData")),
      textInput("out_prefix", "Text you want to prefix output", value = dir_run, width = NULL, placeholder = NULL)
    ),
    mainPanel(),
))


server <- function(input, output) {
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "con_test.csv"
    },
    content = function(file){
      file.copy(file.path(work, "con_test.csv"), file)
    }
  )
}

shinyApp(ui, server)
