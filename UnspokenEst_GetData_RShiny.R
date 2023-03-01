if (exists("env_shiny")) rm(env_shiny)
env_shiny <- new.env(parent = emptyenv())

# Define UI for data upload app ----
env_shiny$ui_unspoken1 <- fluidPage(
  # App title ----
  titlePanel("Upload files for Unspoken into R"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      fileInput("file_attraction", "Attraction Data (csv or rds)",
                multiple = FALSE,
                accept = c(".csv",".rds")),
      fileInput("file_conv","Conversion Data (csv or rds)",
                multiple = FALSE,
                accept = c(".csv",".rds")),
      checkboxInput("use_cov", "Use Covariates (Optional)"),
      conditionalPanel(
        condition = "input.use_cov == true",
        fileInput("file_cov", "Optional: File with Covariates",
                  multiple = FALSE,
                  accept = c(".csv",".rds"))),
      tags$hr(),
      # Checkboxes
      textInput("out_prefix", "Text you want to prefix output", value = "MyUnspoken", width = NULL, placeholder = NULL),
      tags$hr(),
      actionButton("setup_ready","Save Changes to R & Exit Upload", class = "btn-primary")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      h4("Attraction Data (first 10 rows):"),
      h6(tableOutput("file1")),
      h4("Conversion Data (first 10 rows):"),
      h6(tableOutput("file2")),
      h4("Optional Covariate Data (first 10 rows):"),
      h6(tableOutput("file3"))
    )
  )
)


env_shiny$server_unspoken1 <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)
  data1 <- reactive({
    if (is.null(input$file_attraction)) {
      result <- NULL
    } else {
      data_file <- input$file_attraction$datapath
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(data_file, as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(data_file)}
      }  
    } 
    return(result)
  })
  
  data2 <- reactive({
    if (is.null(input$file_conv)) {
      result <- NULL
    } else {
      data_file <- input$file_conv$datapath
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(data_file, as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(data_file)}
      }  
    } 
    return(result)
  })
  
  data3 <- reactive({
    if (is.null(input$file_cov)) {
      result <- NULL
    } else {
      data_file <- input$file_cov$datapath
      ftype <- toupper(substr(data_file,nchar(data_file)-2, nchar(data_file)))
      if (ftype %in% c("CSV", "RDS")){
        if (ftype == "CSV") {result <- read.csv(data_file, as.is=TRUE, check.names = FALSE)}
        if (ftype == "RDS") {result <- readRDS(data_file)}
      }
    }   
    return(result)
  })
  
  output$file1 <- renderTable({data1()[1:10,]})
  output$file2 <- renderTable({data2()[1:10,]})
  output$file3 <- renderTable({data3()[1:10,]})
  
  observeEvent(input$setup_ready, {
    .GlobalEnv$out_prefix <- input$out_prefix
    if (!is.null(data1())) {.GlobalEnv$att_data <- data1()}
    if (!is.null(data2())) {.GlobalEnv$conv_data <- data2()}
    if (input$use_cov){
      .GlobalEnv$data_cov <- data3()
    } else .GlobalEnv$data_cov <- NULL
    stopApp()
  })
}
attach(env_shiny)