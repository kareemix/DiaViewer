library(shiny)
library(bslib)

# Define UI ----
ui <- page_sidebar(
    title = "DIA Viewer",
    selectInput("filename",
        label = h3("Select File:"),
        choices = list.files(path = "./data")
    )
)

# Define server logic ----
server <- function(input, output) {

}

# Run the app ----
shinyApp(ui = ui, server = server)
