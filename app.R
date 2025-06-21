library(shiny)
library(bslib)
library(arrow)
library(ggplot2)

options(shiny.port = 3030)


card_ui <- function(id) {
    ns <- NS(id)
    card(
        layout_columns(
            tags$div(
                selectizeInput(ns("file_name"),
                    label = h5("Select File:"),
                    choices = list.files(path = "./data", pattern = ".*\\.xic\\.parquet$"),
                ),
            ),
            tags$div(
                selectizeInput(
                    ns("peptide_name"),
                    label = h5("Select Peptide:"),
                    choices = NULL,
                    options = list(maxOptions = 1000000000)
                ),
            ),
        ),
        textOutput(ns("loading_text")),
        plotOutput(ns("chart"))
    )
}
card_server <- function(id) {
    moduleServer(id, function(input, output, session) {
        list_peptides <- NULL
        output$chart <- NULL
        get_df <- reactive({
            req(input$file_name)
            read_parquet(
                paste0("./data/", input$file_name),
                col_select = c("pr", "feature", "rt", "value")
            )
        })
        get_pep <- reactive({
            output$loading_text <- renderText("Loaded")
            req(input$peptide_name)
            list_peptides[[input$peptide_name]]
        })
        observeEvent(input$file_name, {
            output$chart <<- NULL
            output$loading_text <- renderText("Loading Parquet...")
            dframe_file <- get_df()
            list_peptides <<- split.data.frame(dframe_file, dframe_file$pr)
            names_peptides <- names(list_peptides)
            updateSelectizeInput(session, "peptide_name", choices = names_peptides, server = TRUE)
            output$loading_text <- renderText("Loading Peptides...")
            observeEvent(input$peptide_name, {
                output$chart <- renderPlot({
                    dframe <- get_pep()
                    dframe_filt <- dframe[dframe$rt != 0, ]
                    rt <- dframe_filt[["rt"]]
                    value <- dframe_filt[["value"]]
                    feature <- dframe_filt[["feature"]]
                    ggplot(data = dframe_filt, mapping = aes(x = rt, y = value, group = feature, color = feature)) +
                        geom_line() +
                        geom_point() +
                        labs(title = paste0("Chromatogram: ", input$peptide_name), x = "Retention Time", y = "Value") +
                        theme_minimal()
                })
            })
        })
    })
}


# Define UI ----
ui <- fluidPage(
    title = "DIA Viewer",
    actionButton("add_card", "Add New Plot"),
    layout_columns(
        div(id = "card_container")
    )
)

# Define server logic ----
server <- function(input, output, session) {
    card_counter <- reactiveVal(0)
    observeEvent(input$add_card, {
        current_card_id <- card_counter() + 1
        card_id <- paste0("item_card_", current_card_id)
        card_counter(current_card_id)
        insertUI(
            selector = "#card_container",
            ui = card_ui(card_id),
            immediate = TRUE
        )
        card_server(card_id)
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
