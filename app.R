library(shiny)
library(bslib)
library(arrow)
library(ggplot2)

options(shiny.port = 3030)


# Define UI ----
ui <- page_sidebar(
    title = "DIA Viewer",
    selectizeInput("file_name",
        label = h5("Select File:"),
        choices = list.files(path = "./data"),
    ),
    selectizeInput(
        "peptide_name",
        label = h5("Select Peptide:"),
        choices = NULL,
        options = list(maxOptions = 1000000000)
    ),
    plotOutput("chart")
)

# Define server logic ----
server <- function(input, output, session) {
    list_peptides <- NULL
    get_df <- reactive({
        req(input$file_name)
        read_parquet(
            paste0("./data/", input$file_name),
            col_select = c("pr", "feature", "rt", "value")
        )
    })
    get_pep <- reactive({
        req(input$peptide_name)
        # print(input$peptide_name)
        list_peptides[[input$peptide_name]]
    })
    observeEvent(get_df(), {
        list_peptides <<- split.data.frame(get_df(), get_df()$pr)
        names_peptides <- names(list_peptides)
        updateSelectizeInput(session, "peptide_name", choices = names_peptides, server = TRUE)
        observeEvent(get_pep(), {
            print(head(get_pep()))
            output$chart <- renderPlot({
                dframe <- get_pep()
                dframe_filt <- dframe[dframe$rt != 0, ]
                ggplot(data = dframe_filt, mapping = aes_string(x = "rt", y = "value")) +
                    geom_line() +
                    labs(title = "Chromatogram", x = "rt", y = "value") +
                    theme_minimal()
            })
        })
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
