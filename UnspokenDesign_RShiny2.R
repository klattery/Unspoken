# user interface
ui <- fluidPage(
  # theme = shinytheme("flatly"),
  titlePanel(title=div(img(src="skim.png", height = 70), "Unspoken Design Generator"),windowTitle = "Unspoken Design"),
  sidebarLayout(
    sidebarPanel(
      h3("Set-Up:"),
      numericInput("numberitems", label = "Select number of total client items:", value = 10),
      numericInput("ntest_perver", label = "Select number of client items to be shown per version (subset):", value = 5),
      numericInput("compitems", label = "Select number of competitor items, always shown in each version, not included in client items:", value = 0),
      h5(strong(htmlOutput("sample_size"))),
      numericInput("show_eachitem_att", label = "Attraction: Number of times to show each item (not a decimal). Recommended is 2+:", value = 2),
      numericInput("show_eachitem", label = "Conversion: Number of times to show each item in conversion (can be a decimal). Recommended is 2.5+:", value = 3),
      numericInput("items_task", label = "Conversion: Number of items to show in each task (not a decimal). 2 is standard:", value = 2),
      numericInput("numberversions", label = "Number of versions:", value = 100),
      fileInput("restrictions", "Choose Restrictions CSV (UTF-8) if there are restrictions on items in versions. CSV should have columns: version, item1, ..., itemn. Must have a row for each version.", multiple = FALSE, accept = c(
        "text/csv",
        "text/comma-separated-values,text/plain",
        ".csv", label = "Restrictions File")),
      fileInput("constraints", "Choose Constraints CSV (UTF-8) if there are items that can't be shown together. CSV should have 2 columns.  Each row has a pair of items that cannot be shown together.", multiple = FALSE, accept = c(
        "text/csv",
        "text/comma-separated-values,text/plain",
        ".csv", label = "Constraints File")),
      actionButton("setup_ready","Setup Done", class = "btn-primary"),
      h4("Download Files once program finishes running:"),
      h3(downloadButton("downconversion","Conversion Design")),
      h3(downloadButton("downattraction","Attraction Design"))
      #h3(downloadButton("downattraction2", "Attraction Design 2"))
      ),
    
    mainPanel(
      wellPanel(style = "background:#E1EBFE",
        h3("Hi, Welcome to the Unspoken-Design Generator Portal!"),
        h4("Please be patient. Your design will take some time to run so we can make sure it is perfectly balanced. 
           Processing time dependes on number of versions. You will know your designs are ready when you see the 1-way tables in the checks below. 
           Have a great day and good luck with your Unspoken Project!"),
        h4("Note: 2 extra concepts will be added to your attraction designs for attention checks.")
        ),
      h4("Conversion Frequency Checks:"),
      h5("Total times each item was shown across all versions:"),
      h6(dataTableOutput("conversion_checks")),
      h5("Total times each item was shown for each version:"),
      h6(dataTableOutput("conversion_checks_2")),
      
      h4("Attraction Frequency Checks:"),
      h5("Total times each item was shown across all versions:"),
      h6(dataTableOutput("att_1_checks")),
      h5("Total times each item was shown for each version:"),
      h6(dataTableOutput("att_1_checks_2")),
    )))

# server
server <- function(input, output) {
  
  # recommended sample size
  output$sample_size <- renderText({
    sample <- ((300*input$numberitems)/input$ntest_perver)
    text<- paste("Your minimum sample size is:","<font color=\"#F27524\">",sample,"</font>", ". Based on your number of total client items and subset size. Sample Size = 300*total_items/items_per_version.")
    text
  })
  

  # read file from user

  restrictions_table <- reactive({
    if (is.null(input$restrictions)) {
      restrictions <- NULL
    } else {
      restrictions <- read.csv(input$restrictions$datapath)
    }
    restrictions})
  
  constraints_table <- reactive({
    if (is.null(input$constraints)) {
      constraints <- NULL
    } else {
      constraints <- read.csv(input$constraints$datapath)
    }
    constraints})
  
  
  observeEvent(input$setup_ready, {
    
    conversion_design <- reactive({
      final <- conversion_function(input$numberitems, input$ntest_perver, input$compitems, input$show_eachitem,  input$items_task, input$numberversions, restrictions_table(), constraints_table())
      final
    })
    
    attraction_design <- reactive({
       final <- attraction_function(conversion_design(), input$show_eachitem_att)
       final
     })
    attraction_design_1 <- reactive({
      final <- attraction_function_1(conversion_design(), input$numberitems, input$compitems, input$ntest_perver, input$show_eachitem_att)
      final
    })
    
    output$conversion_checks <- renderDataTable({dcast(conversion_design(), . ~ item, value.var = 'item', fun.aggregate = length)})
    output$conversion_checks_2 <- renderDataTable({dcast(conversion_design(), version ~ item, value.var = 'item', fun.aggregate = length)})
    
    output$att_1_checks <- renderDataTable({dcast(attraction_design(), . ~ item, value.var = 'item', fun.aggregate = length)})
    output$att_1_checks_2 <- renderDataTable({dcast(attraction_design(), version ~ item, value.var = 'item', fun.aggregate = length)})

    
    output$downconversion <- downloadHandler(
      filename = "ConversionDesign.csv", content = function (file) {write.csv(conversion_design(),file, row.names = FALSE)})
    output$downattraction <- downloadHandler(
       filename = "AttractionDesign.csv", content = function (file) {write.csv(attraction_design(),file, row.names = FALSE)})
  })
}

