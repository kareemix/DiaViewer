library(shiny)
library(bslib)
library(arrow)
library(ggplot2)

options(shiny.port = 3030)


list_peptides <- c()
list_peptides_list <- c()
list_mod_obs <- c()


plot_ui <- function(id, file_text) {
    ns <- NS(id)
    div(
        id = ns("plot_ui"),
        card(
            file_text,
            plotOutput(ns("chart"))
        )
    )
}

plot_server <- function(id, index, selector, yfilter, xfilter) {
    moduleServer(id, function(input, output, session) {
        get_pep <- reactive({
            req(selector())
            return((list_peptides_list[[index]])[[selector()]])
        })
        obs <- observeEvent(selector(), {
            dframe <- get_pep()
            output$chart <- renderPlot({
                dframe_filt <- dframe[dframe$rt != 0, ]
                rt <- dframe_filt[["rt"]]
                value <- dframe_filt[["value"]]
                feature <- dframe_filt[["feature"]]
                pep_plot <- ggplot(data = dframe_filt, mapping = aes(x = rt, y = value, group = feature, color = feature)) +
                    geom_line() +
                    geom_point(size = 3) +
                    labs(title = paste0("Chromatogram: ", selector()), x = "Retention Time", y = "Value") +
                    theme_minimal()
                if (length(dframe_filt$rt) == 0) {
                    pep_plot <- pep_plot + annotate("text",
                        x = 1, y = 1, label = "No Data",
                        color = "black",
                        size = 15,
                    )
                } else {
                    ymax <- max(dframe_filt$value)
                    xmax <- max(dframe_filt$rt)
                    xmin <- dframe_filt$rt[1]
                    xl <- c(xfilter()[1] / 100 * (xmax - xmin) + xmin, xfilter()[2] / 100 * (xmax - xmin) + xmin)
                    yl <- c(yfilter()[1] / 100 * ymax, yfilter()[2] / 100 * ymax)
                    pep_plot <- pep_plot + coord_cartesian(xlim = xl, ylim = yl)
                }
                pep_plot
            })
        })
        obs
    })
}

# Define UI ----
ui <- fluidPage(
    title = "DIA Viewer",
    layout_columns(
        tags$div(
            selectizeInput("file_name",
                label = h5("Select File:"),
                choices = list.files(path = "./data", pattern = ".*\\.xic\\.parquet$"),
                multiple = TRUE,
                selected = list.files(path = "./data", pattern = ".*\\.xic\\.parquet$"),
                # TODO - new card_ui func to add same num of plots as files
            ),
        ),
        tags$div(
            selectizeInput(
                "peptide_name",
                label = h5("Select Peptide:"),
                choices = NULL,
                options = list(maxOptions = 1000000000)
            ),
        ),
    ),
    actionButton("plot_button", label = "Plot"),
    layout_columns(
        sliderInput("yfilter", label = h5("Filter by Value"), min = 0, max = 100, value = c(0, 100)),
        sliderInput("xfilter", label = h5("Filter by Retention Time"), min = 0, max = 100, value = c(0, 100)),
    ),
    div(id = "chart_container")
)

append_unique_list <- function(target, source) {
    unique(append(target, source))
}


# Define server logic ----
server <- function(input, output, session) {
    get_df <- function(name) {
        read_parquet(
            paste0("./data/", name),
            col_select = c("pr", "feature", "rt", "value")
        )
    }
    observeEvent(input$plot_button, {
        removeUI(selector = "#chart_container > *", multiple = TRUE, immediate = TRUE)
        for (i in seq_along(list_mod_obs)) {
            (list_mod_obs[[i]])$destroy()
        }
        updateSelectizeInput(session, "peptide_name", choices = NULL, server = TRUE)
        list_peptides <<- c()
        list_peptides_list <<- c()
        file_list <- c()
        withProgress(message = "Processing Parquet(s)...", {
            for (i in seq_along(input$file_name)) {
                parq <- input$file_name[[i]]
                dframe_file <- get_df(parq)
                file_list <- append(file_list, list(dframe_file))
                incProgress(1 / length(input$file_name))
            }
        })
        names_peptides <- c()
        withProgress(message = "Processing Peptides...", {
            for (i in seq_along(file_list)) {
                dframe_file <- file_list[[i]]
                list_peptides <<- split.data.frame(dframe_file, dframe_file$pr)
                names_peptides <- append_unique_list(names_peptides, names(list_peptides))
                list_peptides_list[length(list_peptides_list) + 1] <<- list(list_peptides)
                incProgress(1 / length(file_list))
            }
        })
        withProgress(message = "Updating Peptide Selection...", {
            updateSelectizeInput(session, "peptide_name", choices = names_peptides, server = TRUE)
        })
        withProgress(message = "Creating Plot(s)...", {
            for (i in seq_along(file_list)) {
                local({
                    this_i <- i
                    new_id <- paste0("plot_", this_i)
                    insertUI(
                        selector = "#chart_container",
                        ui = plot_ui(id = new_id, file_text = input$file_name[this_i]),
                        immediate = TRUE
                    )
                    obs <- plot_server(id = new_id, index = this_i, selector = reactive({
                        input$peptide_name
                    }), yfilter = reactive({
                        input$yfilter
                    }), xfilter = reactive({
                        input$xfilter
                    }))
                    list_mod_obs <<- append(list_mod_obs, list(obs))
                })
                incProgress(1 / length(file_list))
            }
        })
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)
