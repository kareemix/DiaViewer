library(shiny)
library(bslib)
library(arrow)
library(ggplot2)
library(reticulate)


source_python("filter.py")

options(shiny.port = 3030)
options(shiny.maxRequestSize = 100 * 1024^2)


list_peptides <- c()
list_peptides_list <- c()
list_mod_obs <- c()
feature_list <- c()
saved_input <- list(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, 8, 1.5)
report_parq <- NULL


plot_ui <- function(id, file_text) {
    ns <- NS(id)
    div(
        id = ns("plot_ui"),
        div(
            hr(),
            file_text,
            plotOutput(ns("chart")),
        )
    )
}


plot_server <- function(id, file_text, index, selector, feature_sel) {
    moduleServer(id, function(input, output, session) {
        get_pep <- reactive({
            req(selector())
            return((list_peptides_list[[index]])[[selector()]])
        })
        obs <- observeEvent(selector(), {
            dframe <- get_pep()
            output$chart <- renderPlot({
                dframe_filt <- dframe[dframe$rt != 0, ]
                dframe_filt <- dframe_filt[dframe_filt$feature %in% feature_sel(), ]
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
    tags$head(
        tags$style(HTML("
            .selectize-input {
            max-height: 150px;
            overflow-y: auto;
            }
            #filler {
            padding-top: 10px;
            }
            #top_widgets {
            background-color: white;
            border-style: double;
            border-width: 5px;
            z-index: 1020;
            }
            #chart_container {
            padding-top: 275px;
            }
            "))
    ),
    div(id = "filler"),
    fixedPanel(
        id = "top_widgets",
        column(
            3,
            tags$div(
                selectizeInput("file_name",
                    label = h5("Select File:"),
                    choices = list.files(path = "./data", pattern = ".*\\.xic\\.parquet$"),
                    multiple = TRUE,
                    selected = list.files(path = "./data", pattern = ".*\\.xic\\.parquet$"),
                ),
                actionButton("plot_button", label = "Plot"),
            ),
        ),
        column(
            3,
            tags$div(
                selectizeInput(
                    "peptide_name",
                    label = h5("Select Peptide:"),
                    choices = NULL,
                    options = list(maxOptions = 1000000000)
                ),
            ),
        ),
        column(
            2,
            tags$div(
                selectizeInput(
                    "feature_select",
                    label = h5("Select Features:"),
                    choices = NULL,
                    multiple = TRUE
                ),
            ),
        ),
        column(
            4,
            fileInput("report_upload", label = h5("report.parquet:"), accept = ".parquet"),
            actionButton("filter_settings", label = "Filter Settings"),
        ),
    ),
    div(id = "chart_container")
)


append_unique_list <- function(target, source) {
    unique(append(target, source))
}


# Define server logic ----
server <- function(input, output, session) {
    updateCheckboxGroupInput(session, "feature_select", selected = feature_list)
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
        withProgress(message = "Processing Parquet(s)", {
            for (i in seq_along(input$file_name)) {
                parq <- input$file_name[[i]]
                dframe_file <- get_df(parq)
                file_list <- append(file_list, list(dframe_file))
                incProgress(1 / length(input$file_name))
            }
        })
        names_peptides <- c()
        withProgress(message = "Processing Peptides", detail = "This may take a while...", {
            for (i in seq_along(file_list)) {
                dframe_file <- file_list[[i]]
                list_peptides <<- split.data.frame(dframe_file, dframe_file$pr)
                names_peptides <- append_unique_list(names_peptides, names(list_peptides))
                list_peptides_list[length(list_peptides_list) + 1] <<- list(list_peptides)
                feature_list <<- append_unique_list(feature_list, names(split.data.frame(dframe_file, dframe_file$feature)))
                incProgress(1 / length(file_list))
            }
        })
        names_peptides <- sort(names_peptides)
        feature_list <- sort(feature_list)
        withProgress(message = "Updating Peptide Selection", {
            updateSelectizeInput(session, "peptide_name", choices = names_peptides, server = TRUE)
            updateSelectizeInput(session, "feature_select", choices = feature_list, server = TRUE, selected = feature_list)
        })
        withProgress(message = "Creating Plot(s)", {
            for (i in seq_along(file_list)) {
                local({
                    this_i <- i
                    new_id <- paste0("plot_", this_i)
                    insertUI(
                        selector = "#chart_container",
                        ui = plot_ui(id = new_id, file_text = input$file_name[this_i]),
                        immediate = TRUE
                    )
                    obs <- plot_server(
                        id = new_id,
                        file_text = input$file_name[this_i],
                        index = this_i,
                        selector = reactive({
                            input$peptide_name
                        }),
                        feature_sel = reactive({
                            input$feature_select
                        })
                    )
                    list_mod_obs <<- append(list_mod_obs, list(obs))
                })
                incProgress(1 / length(file_list))
            }
        })
    })
    observeEvent(input$filter_settings, {
        showModal(
            modalDialog(
                title = "Filter Settings",
                footer = actionButton("filter_close", label = "Close"),
                fade = FALSE,
                checkboxInput("filter_1", label = "empirical_lib:", value = TRUE),
                checkboxInput("filter_2", label = "peptidoform_mode:", value = TRUE),
                checkboxInput("filter_3", label = "plexdia:", value = FALSE),
                checkboxInput("filter_4", label = "PGMaxLFQ:", value = FALSE),
                checkboxInput("filter_5", label = "QQ:", value = FALSE),
                checkboxInput("filter_6", label = "avg_quality_filter:", value = FALSE),
                checkboxInput("filter_7", label = "filter_peak_width:", value = FALSE),
                numericInput("filter_8", label = h5("min_points_across_peak"), value = 8),
                numericInput("filter_9", label = h5("duty_cycle"), value = 1.5),
            )
        )
    })
    observeEvent(input$filter_close, {
        saved_input <<- list(
            input$filter_1,
            input$filter_2,
            input$filter_3,
            input$filter_4,
            input$filter_5,
            input$filter_6,
            input$filter_7,
            input$filter_8,
            input$filter_9
        )
        if(!is.null(input$report_upload)) {
            report_parq <<- filter_diann(
                report_parq,
                empirical_lib = saved_input[[1]],
                peptidoform_mode = saved_input[[2]],
                plexdia = saved_input[[3]],
                PGMaxLFQ = saved_input[[4]],
                QQ = saved_input[[5]],
                avg_quality_filter = saved_input[[6]],
                filter_peak_width = saved_input[[7]],
                min_points_across_peak = saved_input[[8]],
                duty_cycle = saved_input[[9]]
            )
        }
        removeModal()
    })
    observeEvent(input$report_upload, {
        report_parq <<- read_parquet(input$report_upload$datapath)
        report_parq <<- filter_diann(
            report_parq,
            empirical_lib = saved_input[[1]],
            peptidoform_mode = saved_input[[2]],
            plexdia = saved_input[[3]],
            PGMaxLFQ = saved_input[[4]],
            QQ = saved_input[[5]],
            avg_quality_filter = saved_input[[6]],
            filter_peak_width = saved_input[[7]],
            min_points_across_peak = saved_input[[8]],
            duty_cycle = saved_input[[9]]
        )
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)
