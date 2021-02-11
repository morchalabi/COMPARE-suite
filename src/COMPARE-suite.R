
require(shiny)
require(shinythemes)
require(shinyjs)
require(shinyFiles)
require(shinyBS)
require(ggExtra)
require(DT)
require(visNetwork)
require(igraph)
require(circlize)
require(inlmisc)
require(writexl)

### load the COMPARE workflow functions into global environment
source("functions/1_import.R")
source("functions/2_spillover_compensation.R")
source("functions/3_signal_drift_correction.R")
source("functions/4_viability_correction.R")
source("functions/5_similarity_matrix_generator.R")
source("functions/6_similarity_matrix_heatmap.R")
source("functions/7_negative_control_outlier_detector.R")
source("functions/8_clustering.R")


########################################################################################################################
### --- module UI functions --- ########################################################################################
########################################################################################################################

### UI for Plots view, calls UI functions for all other visualization tabs in Plots view
# @param id String, id of module to be passed to corresponding Server function
plots_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    tabsetPanel(
      
      tabPanel(
        title = "UMAP",
        br(),
        
        Vis_umap_UI(ns("Vis_umap_UI_1"))
        
      ), # end tabPanel
      
      tabPanel(
        title = "Dispersion network",
        br(),
        
        Vis_disp_network_UI(ns("Vis_disp_network_UI_1"))
        
      ), # end tabPanel
      
      tabPanel(
        title = "Samples network",
        br(),
        
        Vis_samples_network_UI(ns("Vis_samples_network_UI_1"))
        
      ) # end tabPanel
      
    ), # end tabsetPanel
    
    fluidRow(
      
      column(
        5,
        
        tabsetPanel(
          
          tabPanel(
            title = "Selected wells",
            br(),
            
            wellPanel(
              style = "border:#f4f4f4",
              
              fluidRow(
                
                column(
                  4,
                  div(actionButton(ns("scatter_wells_select_all_rows"), "Select all", width = "100%"), style = "margin-bottom:20px")
                ), # end column
                
                column(
                  4,
                  div(actionButton(ns("scatter_wells_clear_selected_rows"), "Clear selection", width = "100%"), style = "margin-bottom:20px")
                ), # end column
                
                column(
                  4,
                  div(actionButton(ns("scatter_wells_remove_all"), "Clear table", width = "100%"), style = "margin-bottom:20px")
                ) # end column
                
              ), # end fluidRow
              
              div(
                style = "margin-bottom: 22px; font-size: 89%; height: 500px; overflow: auto;",
                dataTableOutput(ns("scatter_wells_dt"))
              ) # end div
              
            ) # end wellPanel
            
          ), # end tabPanel
          
          tabPanel(
            title = "Control cluster wells",
            br(),
            
            wellPanel(
              style = "border:#f4f4f4",
              
              fluidRow(
                column(
                  4,
                  div(actionButton(ns("ctr_wells_select_all_rows"), "Select all", width = "100%"), style = "margin-bottom:20px")
                ), # end column
                column(
                  4,
                  div(actionButton(ns("ctr_wells_clear_selected_rows"), "Clear selection", width = "100%"), style = "margin-bottom:20px")
                ) # end column
              ), # end fluidRow
              
              div(
                style = "margin-bottom: 22px; font-size: 89%; height: 500px; overflow: auto;",
                dataTableOutput(ns("ctr_wells_dt"))
              ) # end div
              
            ) # end wellPanel
            
          ) # end tabPanel
          
        ), # end tabsetPanel
        
        wellPanel(
          style = "border:#f4f4f4",
          
          fluidRow(
            
            column(
              6,
              div(selectInput(ns("annot_sampling_variable"), "Variable:",
                              choices = c("community","file","drug","control"),
                              width = "100%"),
                  style = "font-size: 95%")
            ), # end column
            
            column(
              6,
              div(numericInput(ns("max_per_group"), "Events per group:", min = 1000, max = 50000, step = 1000, value = 5000),
                  style = "font-size:95%")
            ) # end column
            
          ), # end fluidRow
          
          fluidRow(
            
            column(
              6,
              div(selectInput(ns("transformation_method"), "Transformation:",
                              choices = c("Logicle" = "logicle", "Log2" = "log2", "Log10" = "log10", "None" = "none"),
                              selected = "logicle", width = "100%"),
                  style = "font-size:95%"),
              bsPopover(ns("transformation_method"), title = NULL,
                        content = "Logicle = linear-like scale for low values and log-like scale for high values. Log2 and log10 = shifted log transformations, nagative values are set to 0 and the log is taken from x + 1. None = no transformation.",
                        placement = "right", trigger = "focus", options = list(container = "body"))
            ), # end column
            
            column(
              2,
              div(checkboxInput(ns("center"), "Center", value = FALSE), style = "font-size:95%; margin-top: 30px"),
              bsPopover(ns("center"), title = NULL,
                        content = "Center with respect to mean.",
                        placement = "bottom", trigger = "focus", options = list(container = "body"))
            ), # end column
            
            column(
              2,
              div(checkboxInput(ns("scale"), "Scale", value = FALSE), style = "font-size:95%; margin-top: 30px"),
              bsPopover(ns("scale"), title = NULL,
                        content = "If centered, scale by standard deviation or by root mean square otherwise.",
                        placement = "bottom", trigger = "focus", options = list(container = "body"))
            ) # end column
            
          ), # end fluidRow
          
          hr(),
          
          fluidRow(
            
            column(
              4,
              actionButton(ns("plot_scatter"), "Plot", width = "100%")
            ) # end column
            
          ), # end fluidRow
          
          div(
            style = "font-size:12px; background: #ffffff; border: 1px;
                         border-style: solid; border-color: #dddddd; border-radius: 4px;
                         padding: 10px; height: 80px; overflow: auto; margin-top: 20px;",
            textOutput(ns("scatter_message_center"))
          ) # end div
          
        ) # end wellPanel
        
      ), # end columns
      
      column(
        7,
        
        tabsetPanel(
          
          tabPanel(
            title = "Scatterplot",
            br(),
            
            fluidRow(
              
              column(
                6,
                div(selectInput(ns("scatter_x"), "Channel X:", choices = NULL, selected = NULL, width = "100%"),
                    style = "font-size:95%;")
              ), # end column
              
              column(
                6,
                div(selectInput(ns("scatter_y"), "Channel Y:", choices = NULL, selected = NULL, width = "100%"),
                    style = "font-size:95%")
              ) # end column
              
            ), # end fluidRow
            
            Vis_scatterplot_UI(ns("Vis_scatterplot_UI_1"))
            
          ), # end tabPanel
          
          tabPanel(
            title = "Violins",
            br(),
            
            Vis_violins_UI(ns("Vis_violins_UI_1"))
            
          ) # end tabPanel
          
        ) # end tabsetPanel
        
      ) # end column
      
    ) # end fluidRow
    
  ) # end tagList
  
} # end plots_UI


### UI for UMAP tab
# @param id String, id of module to be passed to corresponding Server function
Vis_umap_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      
      column(
        4,
        div(selectInput(ns("umap_channel"), label = "Channel:", choices = NULL, selected = NULL, width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        2,
        div(selectInput(ns("umap_x"), "X axis:", choices = NULL, selected = NULL, width = "100%"), style = "font-size:95%")
      ),
      
      column(
        2,
        div(selectInput(ns("umap_y"), "Y axis:", choices = NULL, selected = NULL, width = "100%"), style = "font-size:95%")
      ), # end column
      
      column(
        2,
        div(selectInput(ns("umap_dot_size"), "Clique size:", choices = c("Events" = "total_live_cells", "Wells" = "length", "None" = "none"), selected = "total_live_cells", width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        2,
        div(selectInput(ns("umap_labels"), "Clique label:", choices = c("Community", "Clique ID", "None"), selected = "Community"),
            style = "font-size:95%")
      ) # end column
      
    ), # end fluidRow
    
    fluidRow(
      column(
        1,
        textOutput(ns("umap_plot_hover_name"))
      ), # end column
      column(
        1,
        textOutput(ns("umap_plot_hover_centr"))
      ), # end column
      column(
        10,
        textOutput(ns("umap_plot_hover_comm"))
      ) # end column
    ), # end fluidRow
    
    textOutput(ns("umap_plot_hover_wells")),
    textOutput(ns("umap_plot_hover_drugs")),
    
    plotOutput(ns("umap_plot"),
               height = "660px",
               click = ns("umap_plot_click"),
               hover = hoverOpts(id = ns("umap_plot_hover"), delay = 15))
    
  ) # end tagList
  
} # end Vis_umap_UI


### UI for Dispersion network tab
# @param id String, id of module to be passed to corresponding Server function
Vis_disp_network_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    tabsetPanel(
      
      tabPanel(
        title = "Network settings",
        br(),
        
        fluidRow(
          
          column(
            2,
            div(selectInput(ns("disp_network_layout"), "Network layout:",
                            choices = c("FR" = "layout_with_fr", "Hierarchical" = "hierarchical"),
                            selected = "layout_with_fr", width = "100%"),
                style = "font-size:95%"),
            bsPopover(ns("disp_network_layout"), title = NULL, content = "FR = Fruchterman-Reingold; Hierarchical = Largest hub node at the top.",
                      placement = "right", trigger = "focus", options = list(container = "body"))
          ), # end column
          
          column(
            2,
            div(selectInput(ns("disp_node_color_metric"), "Node color:",
                            choices = c("MFI" = "mfi", "Degree" = "degree", "Cells" = "cells",
                                        "Sim vs ctr" = "sim_vs_ctr", "Sim vs all" = "sim_vs_all", "None" = "none"),
                            selected = "mfi"),
                style = "font-size:95%")
          ), # end column
          
          column(
            3,
            div(selectInput(ns("disp_node_color_channel"), label = "Channel for node MFI:",
                            choices = NULL, selected = NULL, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            2,
            div(selectInput(ns("disp_node_label"), "Node label:",
                            choices = c("Well" = "well", "Drug" = "drug", "Community" = "comm"),
                            selected = "well"),
                style = "font-size:95%")
          ), # end column
          
          column(
            2,
            div(selectInput(ns("disp_edge_color_metric"), "Edge color:",
                            choices = c("None" = "none", "Weight" = "weight"),
                            selected = "none"),
                style = "font-size:95%")
          ), # end column
          
          column(
            1,
            div(checkboxInput(ns("disp_highlight_nearest"), "Highlight nearest neighbors", value = FALSE, width = "100%"),
                style = "font-size:95%; margin-top:-4px")
          ), # end column
          
        ) # end fluidRow
        
      ), # end tabPanel
      
      tabPanel(
        title = "Select nodes",
        br(),
        
        fluidRow(
          
          column(
            4,
            div(selectInput(ns("disp_select_by_id"), label = "Node ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            4,
            div(selectInput(ns("disp_select_by_drug"), label = "Drug:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            4,
            div(selectInput(ns("disp_select_by_comm"), label = "Community:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ) # end column
          
        ) # end fluidRow
        
      ) # end tabPanel
      
    ), # end tabsetPanel
    
    fluidRow(
      
      column(
        2,
        div(actionButton(ns("disp_network_add_selected"), "Add to selected wells" , width = "100%"), style = "margin-bottom: 20px")
      ), # end column
      
      column(
        10,
        div(
          style = "background: #ffffff; margin-bottom: 20px; padding: 0px;
                   border: 1px; border-style: solid; border-color: #cccccc; border-radius: 4px;",
          visNetworkOutput(ns("disp_network"), height = 660)
        ) # end div
      ) # end column
      
    ) # end fluidRow
    
  ) # end tagList
  
} # end Vis_disp_network_UI


### UI for Samples network tab
# @param id String, id of module to be passed to corresponding Server function
Vis_samples_network_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    tabsetPanel(
      
      tabPanel(
        title = "Network settings",
        br(),
        
        fluidRow(
          
          column(
            2,
            div(selectInput(ns("samples_network_layout"), "Network layout:",
                            choices = c("FR" = "layout_with_fr", "KK" = "layout_with_kk"),
                            selected = "layout_with_fr"),
                style = "font-size:95%"),
            bsPopover(ns("samples_network_layout"), title = NULL, content = "FR = Fruchterman-Reingold; KK = Kamada-Kawai",
                      placement = "right", trigger = "hover", options = list(container = "body"))
          ), # end column
          
          column(
            2,
            div(selectInput(ns("samples_node_color_metric"), "Node color:",
                            choices = c("MFI" = "mfi", "Degree" = "degree", "Cells" = "cells",
                                        "Sim vs ctr" = "sim_vs_ctr", "Sim vs all" = "sim_vs_all", "None" = "none"),
                            selected = "mfi"),
                style = "font-size:95%")
          ), # end column
          
          column(
            3,
            div(selectInput(ns("samples_node_color_channel"), "Channel for node MFI:", choices = NULL, multiple = FALSE, selected = NULL, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            2,
            div(selectInput(ns("samples_node_label"), "Node label:", choices = c("Well" = "well", "Drug" = "drug", "Community" = "comm"), selected = "well"),
                style = "font-size:95%")
          ), # end column
          
          column(
            2,
            div(selectInput(ns("samples_edge_color_metric"), "Edge color:", choices = c("None" = "none", "Weight" = "weight"), selected = "none"),
                style = "font-size:95%")
          ), # end column
          
          column(
            1,
            div(checkboxInput(ns("samples_highlight_nearest"), "Highlight nearest neighbors", value = FALSE, width = "100%"),
                style = "font-size:95%; margin-top:-4px")
          ), # end column
          
        ) # end fluidRow
        
      ), # end tabPanel
      
      tabPanel(
        title = "Select nodes",
        br(),
        
        fluidRow(
          
          column(
            3,
            div(selectInput(ns("samples_select_by_id"), label = "Node ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            3,
            div(selectInput(ns("samples_select_by_drug"), label = "Drug:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            3,
            div(selectInput(ns("samples_select_by_comm"), label = "Community:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            3,
            div(selectInput(ns("samples_select_by_clique"), label = "Clique ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                style = "font-size:95%")
          ) # end column
          
        ) # end fluidRow
        
      ), # end tabPanel
      
      tabPanel(
        title = "Filter nodes",
        br(),
        
        fluidRow(
          
          column(
            4,
            div(sliderInput(ns("node_degree_range"), "Node degree:", min = 1, max = 100, value = c(1,100), step = 1, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            4,
            div(sliderInput(ns("sim_vs_ctr_range"), "Similarity vs control:", min = 0, max = 100, value = c(0,100), step = 0.1, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            4,
            div(sliderInput(ns("sim_vs_all_range"), "Similarity vs all:", min = 0, max = 100, value = c(0,100), step = 0.1, width = "100%"),
                style = "font-size:95%")
          ) # end column
          
        ) # end fluidRow
        
      ) # end tabPanel
      
    ), # end tabsetPanel
    
    fluidRow(
      
      column(
        2,
        div(actionButton(ns("samples_network_add_selected"), "Add to selected wells" , width = "100%"), style = "margin-bottom: 20px")
      ), # end column
      
      column(
        10,
        div(
          style = "background: #ffffff; margin-bottom: 20px; padding: 0px;
                   border: 1px; border-style: solid; border-color: #cccccc; border-radius: 4px;",
          visNetworkOutput(ns("samples_network"), height = 660)
        ) # end div
      ) # end column
      
    ) # end fluidRow
    
  ) # end tagList
  
} # end Vis_samples_network_UI


### UI for Scatterplot tab
# @param id String, id of module to be passed to corresponding Server function
Vis_scatterplot_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      
      column(
        3,
        div(selectInput(ns("color"), "Point color:", choices = c("Group" = "group", "Density" = "density", "None" = "none"), selected = "group", width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        3,
        div(selectInput(ns("margins"), "Margins:", choices = c("Density" = "density", "Boxplot" = "boxplot", "None" = "none"),
                        selected = "density", width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        3,
        div(selectInput(ns("height"), "Height:", choices = seq(300, 1500, 50), selected = 800, width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        3,
        div(selectInput(ns("point_size"), "Point size:", choices = seq(0.1, 2, 0.1), selected = 0.5, width = "100%"),
            style = "font-size:95%")
      ) # end column
      
    ), # end fluidRow
    
    plotOutput(ns("scatterplot"), height = "700px")
    
  ) # end tagList
  
} # end Vis_scatterplot_UI


### UI for Violins tab
# @param id String, id of module to be passed to corresponding Server function
Vis_violins_UI <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    fluidRow(
      
      column(
        3,
        div(selectInput(inputId = ns("color"), "Violins color:", choices = c("Group" = "group", "None" = "none"),
                        selected = "group", width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        7,
        div(selectInput(ns("channels"), "Channels:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
            style = "font-size:95%")
      ), # end column
      
      column(
        2,
        div(selectInput(ns("width_type"), "Width:", choices = c("Auto" = "auto", "Wide" = "wide"), selected = "auto", width = "100%"),
            style = "font-size:95%")
      ) # end column
      
    ), # end fluidRow
    
    div(
      style = "height: 880px; overflow-x: auto;",
      plotOutput(ns("violins_plot"), height = "860px")
    ) # end div
    
  ) # end tagList
  
} # end Vis_violins_UI


########################################################################################################################
### --- module server functions --- ####################################################################################
########################################################################################################################

### Server for Plots view
# @param id String, id of corresponding plots_UI module
# @param plots_set Reactive expression containing a list of data frames for visualizations in Plots view
# @param paths Reactive expression containing a list of paths to data and out folders
plots_S <- function(id, plots_set, paths) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the plots_set
      plots_set_ref <- reactive({
        validate(need(plots_set(), message = FALSE))
        return(plots_set())
      }) # end plots_set_ref
      
      
      rvs <- reactiveValues(
        clicked_wells = NULL
      ) # end rvs
      
      
      ### update UI controls when a new plots_set is received
      observe({
        
        active_channels <- plots_set_ref()$channels_table[which(plots_set_ref()$channels_table$processed_for_plots == "yes"),]
        
        channel_choices <- active_channels$name
        names(channel_choices) <- active_channels$desc
        
        updateSelectInput(session, "scatter_x", choices = channel_choices)
        updateSelectInput(session, "scatter_y", choices = channel_choices, selected = channel_choices[2])
        
      }) # end observe
      
      
      ### select or deselect rows in the control wells table
      scatter_wells_dt_proxy <- dataTableProxy("scatter_wells_dt")
      
      observeEvent(input$scatter_wells_select_all_rows, {
        selectRows(scatter_wells_dt_proxy, c(input$scatter_wells_dt_rows_all, input$scatter_wells_dt_rows_selected))
      }) # end observeEvent
      
      observeEvent(input$scatter_wells_clear_selected_rows, {
        selectRows(scatter_wells_dt_proxy, NULL)
      }) # end observeEvent
      
      ### clear selected wells
      observeEvent(input$scatter_wells_remove_all, {
        rvs$clicked_wells <- NULL
      }) # end observeEvent
      
      
      ### render the drugs table
      output$ctr_wells_dt <- renderDT(
        plots_set_ref()$ctr_wells[,c("file","drug","concentration","control","live_cells")],
        filter = "top",
        rownames = FALSE,
        options = list(pageLength = 500, lengthMenu = c(10,25,50,100,500,1000), searchHighlight = TRUE)
      ) # end output$ctr_wells_dt
      
      
      ### select or deselect rows in the control wells table
      ctr_wells_dt_proxy <- dataTableProxy("ctr_wells_dt")
      
      observeEvent(input$ctr_wells_select_all_rows, {
        selectRows(ctr_wells_dt_proxy, c(input$ctr_wells_dt_rows_all, input$ctr_wells_dt_rows_selected))
      }) # end observeEvent
      
      observeEvent(input$ctr_wells_clear_selected_rows, {
        selectRows(ctr_wells_dt_proxy, NULL)
      }) # end observeEvent
      
      
      ### UMAP module
      umap_clicked <- Vis_umap_S("Vis_umap_UI_1", plots_set_ref)
      
      
      ### Dispersion network module
      disp_network_clicked <- Vis_disp_network_S("Vis_disp_network_UI_1", plots_set_ref)
      
      
      ### Samples network module
      samples_network_clicked <- Vis_samples_network_S("Vis_samples_network_UI_1", plots_set_ref)
      
      
      ### collect clicked wells from all visualization modules
      observe({
        rvs$clicked_wells <- union(isolate(rvs$clicked_wells),
                                   c(umap_clicked()$wells,
                                     disp_network_clicked()$wells,
                                     samples_network_clicked()$wells))
      }) # end observe
      
      
      ### subset of annotation with wells clicked in umap and network plots
      scatter_wells <- reactive({
        if (is.null(rvs$clicked_wells)) {
          return(NULL)
        } else {
          rows_to_select <- which(plots_set_ref()$drugs$file %in% rvs$clicked_wells & !(plots_set_ref()$drugs$file %in% plots_set_ref()$ctr_wells$file))
          validate(need(rows_to_select, message = FALSE))
          return(plots_set_ref()$drugs[rows_to_select,])
        } # end else
      }) # end scatter_wells
      
      
      ### observer to render the datatable with wells selected in umap and network plots
      observeEvent(scatter_wells(), {
        
        if (is.null(scatter_wells())) {
          output$scatter_wells_dt <- renderDT(NULL)
        } else {
          output$scatter_wells_dt <- renderDT(
            scatter_wells()[,c("file","drug","concentration","community","live_cells")],
            filter = "top",
            rownames = FALSE,
            options = list(pageLength = 500, lengthMenu = c(10,25,50,100,500,1000), searchHighlight = TRUE)
          ) # end output$scatter_wells_dt
        } # end else
        
      }) # end observe
      
      
      selected_sampling_method <- reactive({
        return(input$annot_sampling_variable)
      }) # end selected_sampling_method
      
      
      sampling_groups <- reactive({
        return(unique(rbind(scatter_wells()[input$scatter_wells_dt_rows_selected,],
                            plots_set_ref()$ctr_wells[input$ctr_wells_dt_rows_selected,])[,selected_sampling_method()]))
      }) # end sampling_groups
      
      
      ### collect data for scatter plot and violins plot modules
      observeEvent(
        {
          input$plot_scatter
          input$scatter_x
          input$scatter_y
        }, {
          
          if (input$scatter_x != input$scatter_y) {
            
            scatter_data <- reactiveVal(event_scatter_data())
            
            # scatterplot module
            Vis_scatterplot_S("Vis_scatterplot_UI_1", scatter_data)
            
          } # end if
          
          violins_data <- reactiveVal(event_violins_data())
          
          # density violins plot module
          Vis_violins_S("Vis_violins_UI_1", violins_data)
          
      }) # end observeEvent
      
      
      ### data for event scatterplot
      event_scatter_data <- reactive({
        
        active_channels <- plots_set_ref()$channels_table[which(plots_set_ref()$channels_table$processed_for_plots == "yes"),]
        
        channel_choices <- active_channels$name
        names(channel_choices) <- active_channels$desc
        
        return(
          list(
            data = sampled_scatterplot_data(),
            selected_var = selected_sampling_method(),
            channel_choices = channel_choices,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y
          ) # end list
        ) # end return
        
      }) # end event_scatter_data
      
      
      ### data for event violins plot
      event_violins_data <- reactive({
        
        active_channels <- plots_set_ref()$channels_table[which(plots_set_ref()$channels_table$processed_for_plots == "yes"),]
        
        channel_choices <- active_channels$name
        names(channel_choices) <- active_channels$desc
        
        return(
          list(
            data = sampled_event_exprs_data(),
            selected_var = selected_sampling_method(),
            channel_choices = channel_choices,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y
          ) # end list
        ) # end return
        
      }) # end event_plots_data
      
      
      ### event expression data after downsampling, filtered for x and y channels for scatterplot
      sampled_scatterplot_data <- reactive({
        
        validate(need(sampled_event_exprs_data(), message = FALSE))
        validate(need(selected_sampling_method(), message = FALSE))
        
        # make copy with selected x and y channels
        d_scatter <- sampled_event_exprs_data()[,c(input$scatter_x, input$scatter_y)]
        
        # attach variable and variable color
        d_scatter$var_name  <- sampled_event_exprs_data()[,selected_sampling_method()]
        d_scatter$var_color <- sampled_event_exprs_data()[,paste0(selected_sampling_method(), "_color")]
        
        # get density in 2d for the x and y points
        d_scatter$density <- get_density(d_scatter[,1], d_scatter[,2], h = c(1,1), n = 100)
        
        html(id = "scatter_message_center", html = paste0(Sys.time(), " : Rendering ...</br>"), add = FALSE)
        
        return(d_scatter)
        
      }) # end sampled_scatterplot_data
      
      
      ### sampled event expression data
      sampled_event_exprs_data <- reactive({
        
        validate(need(event_exprs_data(), message = FALSE))
        validate(need(selected_sampling_method(), message = FALSE))
        
        html(id = "scatter_message_center", html = paste0(Sys.time(), " : Sampling ...</br>"), add = FALSE)
        
        d <- event_exprs_data()[which(event_exprs_data()[,selected_sampling_method()] %in% sampling_groups()),]
        
        sv <- selected_sampling_method()
        
        var_groups <- unique(d[,sv])
        
        d_sampled <- NULL
        
        for (group in var_groups) {
          
          d_group <- d[which(d[,sv] == group),]
          
          if (input$max_per_group >= nrow(d_group)) {
            d_sampled <- rbind(d_sampled, d_group)
          } else {
            d_group_sampled <- d_group[sample(nrow(d_group), input$max_per_group),]
            d_group_sampled <- d_group_sampled[order(as.numeric(rownames(d_group_sampled))),]
            d_sampled <- rbind(d_sampled, d_group_sampled)
          } # end else
          
        } # end for
        
        #write.table(d_sampled, "../event_exprs_data_sampled.txt", sep = "\t", quote = FALSE)
        
        return(d_sampled)
        
      }) # end sampled_event_exprs_data
      
      
      ### event expression data
      event_exprs_data <- reactive({
        
        # collect highlighted rows from selected wells and control wells tables
        wells <- rbind(scatter_wells()[input$scatter_wells_dt_rows_selected,], plots_set_ref()$ctr_wells[input$ctr_wells_dt_rows_selected,])
        
        if (nrow(wells) == 0) {
          html(id = "scatter_message_center",
               html = "To generate Scatterplot and Violins, click to add wells from UMAP plot, Dispersion,
                       or Samples networks and highlight rows in Selected wells or Control cluster wells.",
               add = FALSE)
          return(NULL)
        } # end if
        
        active_channels <- plots_set_ref()$channels_table[which(plots_set_ref()$channels_table$processed_for_plots == "yes"),]
        
        channel_names <- active_channels$name
        
        files <- paste0(paths()$data, wells$file, ".fcs")
        
        # read fcs files into FlowSet
        html(id = "scatter_message_center", html = paste0(Sys.time(), " : Reading ", length(files), " FCS files ...</br>"), add = FALSE)
        fset_obj <- read.flowSet(files, transformation = FALSE)
        
        # logicle transformation if selected
        if (input$transformation_method == "logicle") {
          html(id = "scatter_message_center", html = paste0(Sys.time(), " : Applying ", input$transformation_method, " transformation ...</br>"), add = TRUE)
          fset_obj <- transform(fset_obj, transformList(channel_names, logicleTransform(w = 0.5, t = 262144, m = 4.5)))
        } # end if
        
        
        vars <- c("community","file","drug","control")
        
        # extract into data frame and record the start and end row positions of each well in the exprs_data
        html(id = "scatter_message_center", html = paste0(Sys.time(), " : Collecting expression data for scatterplot ...</br>"), add = TRUE)
        
        exprs_data <- NULL
        var_cols <- NULL
        wells$start <- NA
        wells$end <- NA
        
        for (w in 1:nrow(wells)) {
          
          # start position
          if (w == 1)
            wells[w,]$start <- 1
          else
            wells[w,]$start <- nrow(exprs_data) + 1
          
          exprs_data <- rbind(exprs_data, fset_obj[[w]]@exprs[,channel_names, drop = FALSE])
          
          # end position
          wells[w,]$end <- nrow(exprs_data)
          
          # columns with specifications of selected variables for each cell
          var_cols_per_well <- NULL
          for (v in vars) {
            var_col <- rep(wells[w,v], wells[w,]$end - wells[w,]$start + 1)
            var_cols_per_well <- cbind(var_cols_per_well, var_col)
          } # end for
          
          var_cols <- rbind(var_cols, var_cols_per_well)
          
        } # end for
        
        
        html(id = "scatter_message_center", html = paste0(Sys.time(), " : Making color palettes for variables ...</br>"), add = TRUE)
        
        # color columns for variables
        var_color_cols <- NULL
        for (c in 1:ncol(var_cols)) {
          color_palette <- GetColors(length(unique(var_cols[,c])), scheme = "smooth rainbow", start = 0.2, end = 0.9, bias = 1.6)
          var_color_cols <- cbind(var_color_cols, color_palette[as.factor(var_cols[,c])])
        } # end for
        
        
        # log transformation if selected
        if (input$transformation_method %in% c("log2","log10")) {
          
          html(id = "scatter_message_center",
               html = paste0(Sys.time(), " : Applying ", input$transformation_method, " transformation, negative values are set to 0 ...</br>"),
               add = TRUE)
          
          if (input$transformation_method == "log2") 
            b = 2
          else if (input$transformation_method == "log10")
            b = 10
          
          # set negative values to 0
          exprs_data[exprs_data < 0] <- 0
          
          # shifted log-transform with base b
          exprs_data <- log(exprs_data + 1, b)
          
        } # end if
        
        # center and scale if selected
        if (input$center & input$scale)
          html(id = "scatter_message_center", html = paste0(Sys.time(), " : Centering and scaling ... </br>"), add = TRUE)
        else if (input$center & !(input$scale))
          html(id = "scatter_message_center", html = paste0(Sys.time(), " : Centering ... </br>"), add = TRUE)
        else if (!(input$center) & input$scale)
          html(id = "scatter_message_center", html = paste0(Sys.time(), " : Scaling ... </br>"), add = TRUE)
        
        exprs_data  <- scale(exprs_data, center = input$center, scale = input$scale)
        
        
        ### compile the final expression data with variables and variable colors
        exprs_data <- as.data.frame(cbind(exprs_data, var_cols, var_color_cols))
        colnames(exprs_data) <- c(channel_names, vars, paste0(vars, "_color"))
        
        for (i in 1:length(channel_names)) {
          exprs_data[,i] <- as.double(exprs_data[,i])
        } # end for
        
        #write.table(exprs_data, "../event_exprs_data.txt", sep = "\t", quote = FALSE)
        
        return(exprs_data)
        
      }) # end event_exprs_data
      
      
    } # end function
    
  ) # end moduleServer
  
} # end plots_S


### Server for UMAP tab
# @param id String, id of corresponding Vis_umap_UI module
# @param plots_set Reactive expression containing a list of data frames for visualizations in Plots view
# @return selected_wells Reactive value containing a list of clicked well identifiers
Vis_umap_S <- function(id, plots_set) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the plots_set
      plots_set_ref <- reactive({
        validate(need(plots_set(), message = FALSE))
        return(plots_set())
      }) # end plots_set_ref
      
      selected_wells <- reactiveVal(NULL)
      
      ### update plot-specific controls that depend on the input plots_set
      observe({
        
        validate(need(plots_set_ref(), message = FALSE))
        
        updateSelectInput(session, "umap_x", choices = colnames(plots_set_ref()$umap))
        updateSelectInput(session, "umap_y", choices = colnames(plots_set_ref()$umap), selected = colnames(plots_set_ref()$umap)[2])
        updateSelectInput(session, "umap_channel", choices = as.character(plots_set_ref()$channels_table[plots_set_ref()$channels,]$desc))
        
      }) # end observe
      
      
      ### placeholders for hover messages
      html(id = "umap_plot_hover_name",
           html = paste("<strong>Clique:</strong>"),
           add = FALSE)
      html(id = "umap_plot_hover_centr",
           html = paste("<strong>MFI:</strong>"),
           add = FALSE)
      html(id = "umap_plot_hover_comm",
           html = paste("<strong>Community:</strong>"),
           add = FALSE)
      html(id = "umap_plot_hover_wells",
           html = paste("<strong>Wells:</strong>"),
           add = FALSE)
      html(id = "umap_plot_hover_drugs",
           html = paste("<strong>Drugs:</strong>"),
           add = FALSE)
      
      
      ### data for the umap scatterplot
      umap_plot_data <- reactive({
        
        validate(need(plots_set_ref(), message = FALSE))
        validate(need(input$umap_x, message = FALSE))
        validate(need(input$umap_y, message = FALSE))
        validate(need(input$umap_channel, message = FALSE))
        
        d <- data.frame(
          x = plots_set_ref()$umap[,input$umap_x],
          y = plots_set_ref()$umap[,input$umap_y],
          centr = plots_set_ref()$centroids[,input$umap_channel],
          name = rownames(plots_set_ref()$umap),
          wells = c("", plots_set_ref()$cliques$File),
          drugs = c("", plots_set_ref()$cliques$Clique),
          comm = c("Control", plots_set_ref()$cliques$Community)
        ) # end data.frame
        
        if (input$umap_dot_size == "total_live_cells") {
          d$size <- c(max(plots_set_ref()$cliques$total_live_cells), plots_set_ref()$cliques$total_live_cells)
          size_title <- "Number\nof events"
        } else if (input$umap_dot_size == "length") {
          d$size <- c(max(plots_set_ref()$cliques$length), plots_set_ref()$cliques$length)
          size_title <- "Number\nof wells"
        } else {
          size_title <- "None"
        } # end else
        
        #write.table(d, "umap_plot_data.txt", sep = "\t", quote = FALSE)
        
        return(list(d = d, x_lab = input$umap_x, y_lab = input$umap_y, channel = input$umap_channel, size_title = size_title))
        
      }) # end umap_plot_data
      
      
      ### render umap scatterplot
      output$umap_plot <- renderPlot({
        
        validate(need(umap_plot_data(), message = FALSE))
        
        d <- umap_plot_data()$d
        
        border_strokes <- rep(0.4, nrow(d))
        border_strokes[which(d$name == "Control")] <- 2
        
        border_cols <- rep("#cccccc", nrow(d))
        border_cols[which(d$name == "Control")] <- "#008000"
        
        node_shapes <- rep(21, nrow(d))
        node_shapes[which(d$name == "Control")] <- 23
        
        if (input$umap_dot_size == "total_live_cells") {
          
          max_size <- max(d$size)
          size_limits <- c(0, max_size)
          size_breaks = seq(0, max_size, max_size / 10)
          size_labels <- c(0, round(seq(max_size / 10, max_size, max_size / 10), 0))
          
          p <- ggplot(d, aes(x = x, y = y)) +
            geom_point(aes(fill = centr, size = size), shape = node_shapes, alpha = 0.8, color = border_cols, stroke = border_strokes) +
            scale_size(breaks = size_breaks, range = c(1,11), labels = size_labels) +
            expand_limits(size = size_limits) +
            labs(x = umap_plot_data()$x_lab, y = umap_plot_data()$y_lab) +
            guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 14)),
                   size = guide_legend(order = 2, title = umap_plot_data()$size_title, label.theme = element_text(size = 14))) +
            theme_light() +
            theme(legend.position = "left", legend.title = element_text(size = 16),
                  axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
            scale_fill_gradientn(colours = c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
                                 limits = c(min(d$centr), max(d$centr)),
                                 breaks = round(seq(min(d$centr) + 0.01, max(d$centr) - 0.01, length.out = 3), 2))
          
        } else if (input$umap_dot_size == "length") {
          
          max_size <- max(d$size)
          size_limits <- c(0, max_size)
          size_breaks = seq(0, max_size, max_size / 5)
          size_labels <- c(0, round(seq(max_size / 5, max_size, max_size / 5), 1))
          
          p <- ggplot(d, aes(x = x, y = y)) +
            geom_point(aes(fill = centr, size = size), shape = node_shapes, alpha = 0.8, color = border_cols, stroke = border_strokes) +
            scale_size(breaks = size_breaks, range = c(1,11), labels = size_labels) +
            expand_limits(size = size_limits) +
            labs(x = umap_plot_data()$x_lab, y = umap_plot_data()$y_lab) +
            guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 14)),
                   size = guide_legend(order = 2, title = umap_plot_data()$size_title, label.theme = element_text(size = 14))) +
            theme_light() +
            theme(legend.position = "left", legend.title = element_text(size = 16),
                  axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
            scale_fill_gradientn(colours = c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
                                 limits = c(min(d$centr), max(d$centr)),
                                 breaks = round(seq(min(d$centr) + 0.01, max(d$centr) - 0.01, length.out = 3), 2))
          
        } else {
          
          p <- ggplot(d, aes(x = x, y = y)) +
            geom_point(aes(fill = centr), size = 10, shape = node_shapes, alpha = 0.8, color = border_cols, stroke = border_strokes) +
            labs(x = umap_plot_data()$x_lab, y = umap_plot_data()$y_lab) +
            guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 16, angle = 45, hjust = 1))) +
            theme_light() +
            theme(legend.position = "left",
                  legend.title = element_text(size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 18)) +
            scale_fill_gradientn(colours = c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
                                 limits = c(min(d$centr), max(d$centr)),
                                 breaks = round(seq(min(d$centr) + 0.01, max(d$centr) - 0.01, length.out = 3), 2))
          
        } # end else
        
        if (input$umap_labels == "Clique ID") {
          
          p <- p + geom_text_repel(box.padding = unit(0.3, "lines"), label = d$name,
                                   size = 5.5, show.legend = FALSE, fontface = 1,
                                   max.overlaps = 50)
          
        } else if (input$umap_labels == "Community") {
          
          color_palette <- GetColors(length(unique(d$comm)), scheme = "jet")
          
          p <- p + geom_text_repel(box.padding = unit(0.6, "lines"), label = d$comm,
                                   size = 7, show.legend = FALSE, fontface = 2, color = color_palette[as.factor(d$comm)],
                                   max.overlaps = 50)
          
        } # end else if
        
        p
        
      }) # end output$umap_plot
      
      
      ### observer to display data on hover point from line plot
      observeEvent(input$umap_plot_hover, {
        
        d <- umap_plot_data()$d
        
        validate(need(d, message = FALSE))
        
        hover <- input$umap_plot_hover
        point <- nearPoints(d, hover, threshold = 10, maxpoints = 1, xvar = "x", yvar = "y")
        
        if (nrow(point) > 0) {
          html(id = "umap_plot_hover_name",
               html = paste("<strong>Clique:</strong>", point$name),
               add = FALSE)
          html(id = "umap_plot_hover_centr",
               html = paste("<strong>MFI:</strong>", round(point$centr, 3)),
               add = FALSE)
          html(id = "umap_plot_hover_comm",
               html = paste("<strong>Community:</strong>", point$comm),
               add = FALSE)
          html(id = "umap_plot_hover_wells",
               html = paste("<strong>Wells:</strong>", point$wells),
               add = FALSE)
          html(id = "umap_plot_hover_drugs",
               html = paste("<strong>Drugs:</strong>", point$drugs),
               add = FALSE)
        } else {
          html(id = "umap_plot_hover_name",
               html = paste("<strong>Clique:</strong>"),
               add = FALSE)
          html(id = "umap_plot_hover_centr",
               html = paste("<strong>MFI:</strong>"),
               add = FALSE)
          html(id = "umap_plot_hover_comm",
               html = paste("<strong>Community:</strong>"),
               add = FALSE)
          html(id = "umap_plot_hover_wells",
               html = paste("<strong>Wells:</strong>"),
               add = FALSE)
          html(id = "umap_plot_hover_drugs",
               html = paste("<strong>Drugs:</strong>"),
               add = FALSE)
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to record clicked wells in plate heatmaps
      observeEvent(input$umap_plot_click, {
        
        d <- umap_plot_data()$d
        
        validate(need(d, message = FALSE))
        
        click <- input$umap_plot_click
        point <- nearPoints(d, click, threshold = 10, maxpoints = 1, xvar = "x", yvar = "y")
        
        if (nrow(point) > 0)
          selected_wells(list(wells = strsplit(point$wells, split = ",")[[1]]))
        
      }) # end observeEvent
      
      return(selected_wells)
      
    } # end function
    
  ) # end moduleServer
  
} # end Vis_umap_S


### Server for Dispersion network tab
# @param id String, id of corresponding Vis_disp_network_UI module
# @param plots_set Reactive expression containing a list of data frames for visualizations in Plots view
# @return selected_wells Reactive value containing a list of clicked well identifiers
Vis_disp_network_S <- function(id, plots_set) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the plots_set
      plots_set_ref <- reactive({
        validate(need(plots_set(), message = FALSE))
        return(plots_set())
      }) # end plots_set_ref
      
      selected_wells <- reactiveVal(NULL)
      
      ### update plot-specific controls that depend on the input plots_set
      observe({
        
        validate(need(plots_set_ref(), message = FALSE))
        
        updateSelectInput(session, "disp_node_color_channel", choices = as.character(plots_set_ref()$channels_table[plots_set_ref()$channels,]$desc))
        
      }) # end observe
      
      ### samples graph igraph object
      disp_igraph <- reactive({
        
        validate(need(plots_set_ref()$drugs, message = FALSE))
        
        nodes <- data.frame(
          id = plots_set_ref()$drugs$file,
          color.border = "#999999",
          borderWidth = 1,
          color.highlight.border = "#d62d20",
          degree = plots_set_ref()$drugs$degree_disp,
          value = plots_set_ref()$drugs$live_cells,
          drug = plots_set_ref()$drugs$drug,
          drug_label = paste(plots_set_ref()$drugs$drug, plots_set_ref()$drugs$concentration, sep = "_"),
          community = plots_set_ref()$drugs$community,
          sim_vs_ctr = plots_set_ref()$drugs$sim_vs_control,
          sim_vs_all = plots_set_ref()$drugs$sim_vs_all
        ) # end data.frame 
        
        # assign node color based on selected metric
        if (input$disp_node_color_metric == "none") {
          
          nodes$color.background <- "#f1f1f1"
          nodes$color.highlight.background <- "#f1f1f1"
          
        } else {
          
          if (input$disp_node_color_metric == "degree")
            col_values <- plots_set_ref()$drugs$degree_disp
          else if (input$disp_node_color_metric == "cells")
            col_values <- plots_set_ref()$drugs$live_cells
          else if (input$disp_node_color_metric == "mfi")
            col_values <- plots_set_ref()$mfis[,input$disp_node_color_channel]
          else if (input$disp_node_color_metric == "sim_vs_ctr")
            col_values <- plots_set_ref()$drugs$sim_vs_control
          else if (input$disp_node_color_metric == "sim_vs_all")
            col_values <- plots_set_ref()$drugs$sim_vs_all
          
          col_fun <- colorRamp2(c(min(col_values), median(col_values), max(col_values)),
                                c("#0057e7", "#fdf498", "#d62d20"))
          nodes$color.background <- col_fun(col_values)
          nodes$color.highlight.background <- col_fun(col_values)
          
        } # end else
        
        # add control node
        nodes <- rbind(nodes, data.frame(id = "Control",
                                         color.border = "#999999",
                                         borderWidth = 1,
                                         color.highlight.border = "#d62d20",
                                         degree = length(which(plots_set_ref()$disp_network$from == "Control")) + length(which(plots_set_ref()$disp_network$to == "Control")),
                                         value = mean(plots_set_ref()$drugs[which(plots_set_ref()$drugs$degree_disp == 0),]$live_cells),
                                         drug = "Control",
                                         drug_label = "Control",
                                         community = 0,
                                         sim_vs_ctr = NA,
                                         sim_vs_all = NA,
                                         color.background = "#f1f1f1",
                                         color.highlight.background = "#f1f1f1"))
        
        # assign node color to the control node based on selected metric
        if (!(input$disp_node_color_metric == "none")) {
          
          if (input$disp_node_color_metric == "degree") {
            
            ctr_value <- nodes[which(nodes$id == "Control"),]$degree
            col_values <- c(plots_set_ref()$drugs$degree_disp, ctr_value)
            
          } else if (input$disp_node_color_metric == "cells") {
            
            ctr_value <- nodes[which(nodes$id == "Control"),]$value
            col_values <- c(plots_set_ref()$drugs$live_cells, ctr_value)
            
          } else if (input$disp_node_color_metric == "mfi") {
            
            ctr_value <- median(plots_set_ref()$mfis[which(plots_set_ref()$drugs$degree_disp == 0),input$disp_node_color_channel])
            col_values <- c(plots_set_ref()$mfis[,input$disp_node_color_channel], ctr_value)
            
          } else if (input$disp_node_color_metric == "sim_vs_ctr") {
            
            ctr_value <- mean(plots_set_ref()$drugs[which(plots_set_ref()$drugs$degree_disp == 0),]$sim_vs_control)
            col_values <- c(plots_set_ref()$drugs$sim_vs_control, ctr_value)
            
          } else if (input$disp_node_color_metric == "sim_vs_all") {
            
            ctr_value <- mean(plots_set_ref()$drugs[which(plots_set_ref()$drugs$degree_disp == 0),]$sim_vs_all)
            col_values <- c(plots_set_ref()$drugs$sim_vs_all, ctr_value)
            
          } # end esle if
          
          col_fun <- colorRamp2(c(min(col_values), median(col_values), max(col_values)),
                                c("#0057e7", "#fdf498", "#d62d20"))
          nodes[which(nodes$id == "Control"),]$color.background <- col_fun(ctr_value)
          nodes[which(nodes$id == "Control"),]$color.highlight.background <- col_fun(ctr_value)
          
        } # end if
        
        
        # assign node label based on user input
        if (input$disp_node_label == "well")
          nodes$label <- nodes$id
        else if (input$disp_node_label == "drug")
          nodes$label <- nodes$drug_label
        else if (input$disp_node_label == "comm")
          nodes$label <- nodes$community
        
        # make edges data frame
        edges <- data.frame(
          from = plots_set_ref()$disp_network$from,
          to = plots_set_ref()$disp_network$to,
          weight = plots_set_ref()$disp_network$weight
        ) # end data.frame
        
        # assign edge colors
        if (input$disp_edge_color_metric == "weight") {
          edges$color <- plots_set_ref()$disp_network$color_weight
        } else {
          edges$color <- "#cbcbcb"
        } # end esle
        
        # filter nodes to only those that are in the edge list
        nodes <- nodes[which(nodes$id %in% unique(c(edges$from, edges$to))),]
        
        # make igraph object
        g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
        
        # update the node id input box
        updateSelectInput(session, "disp_select_by_id", choices = V(g)$name[order(V(g)$name)])
        updateSelectInput(session, "disp_select_by_drug", choices = V(g)$drug[order(V(g)$drug)])
        updateSelectInput(session, "disp_select_by_comm", choices = V(g)$community[order(as.integer(V(g)$community))])
        
        return(g)
        
      }) # end disp_igraph
      
      
      ### render the dispersion network
      output$disp_network <- renderVisNetwork({
        
        g <- disp_igraph()
        
        if (is.null(g))
          return(NULL)
        
        # render the network from igraph object
        if (input$disp_network_layout == "layout_with_fr") {
          
          visIgraph(g, idToLabel = FALSE, layout = "layout_with_fr", randomSeed = 211) %>%
            visOptions(highlightNearest = input$disp_highlight_nearest, autoResize = TRUE) %>%
            visInteraction(multiselect = TRUE)
          
        } else if (input$disp_network_layout == "hierarchical") {
          
          visIgraph(g, idToLabel = FALSE) %>%
            visHierarchicalLayout(direction = "UD", levelSeparation = 130, nodeSpacing = 120, sortMethod = "hubsize") %>%
            visOptions(highlightNearest = input$disp_highlight_nearest, autoResize = TRUE) %>%
            visInteraction(multiselect = TRUE)
          
        } # end else if
        
      }) # end output$disp_network
      
      
      ### observer to select nodes by id
      observeEvent(input$disp_select_by_id, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$disp_select_by_id))
          visNetworkProxy(ns("disp_network")) %>% visUnselectAll()
        else
          visNetworkProxy(ns("disp_network")) %>% visSelectNodes(id = input$disp_select_by_id)
        
      }) # end observeEvent
      
      
      ### observer to select nodes by drug
      observeEvent(input$disp_select_by_drug, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$disp_select_by_drug))
          visNetworkProxy(ns("disp_network")) %>% visUnselectAll()
        else {
          g <- disp_igraph()
          visNetworkProxy(ns("disp_network")) %>%
            visSelectNodes(id = V(g)$name[which(V(g)$drug %in% input$disp_select_by_drug)])
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to select nodes by community
      observeEvent(input$disp_select_by_comm, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$disp_select_by_comm))
          visNetworkProxy(ns("disp_network")) %>% visUnselectAll()
        else {
          g <- disp_igraph()
          visNetworkProxy(ns("disp_network")) %>%
            visSelectNodes(id = V(g)$name[which(V(g)$community %in% input$disp_select_by_comm)]) 
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to add selected nodes form the dispersion network into the scatter wells table
      observe({
        
        input$disp_network_add_selected
        
        ns <- session$ns
        
        visNetworkProxy(ns("disp_network")) %>% visGetSelectedNodes()
        
        validate(need(input$disp_network_selectedNodes, message = FALSE))
        
        selected_wells(list(wells = setdiff(input$disp_network_selectedNodes, "Control")))
        
      }) # end observe
      
      return(selected_wells)
      
    } # end function
    
  ) # end moduleServer
  
} # end Vis_disp_network_S


### Server for Samples network tab
# @param id String, id of corresponding Vis_samples_network_UI module
# @param plots_set Reactive expression containing a list of data frames for visualizations in Plots view
# @return selected_wells Reactive value containing a list of clicked well identifiers
Vis_samples_network_S <- function(id, plots_set) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the plots_set
      plots_set_ref <- reactive({
        validate(need(plots_set(), message = FALSE))
        return(plots_set())
      }) # end plots_set_ref
      
      selected_wells <- reactiveVal(NULL)
      
      ### update plot-specific controls that depend on the input plots_set
      observe({
        
        validate(need(plots_set_ref(), message = FALSE))
        
        updateSelectInput(session, "samples_node_color_channel", choices = as.character(plots_set_ref()$channels_table[plots_set_ref()$channels,]$desc))
        
      }) # end observe
      
      
      ### samples graph igraph object
      samples_igraph <- reactive({
        
        validate(need(plots_set_ref()$drugs, message = FALSE))
        
        nodes <- data.frame(
          id = plots_set_ref()$drugs$file,
          color.border = "#999999",
          borderWidth = 1,
          color.highlight.border = "#d62d20",
          value = plots_set_ref()$drugs$live_cells,
          degree = plots_set_ref()$drugs$degree_samples,
          sim_vs_ctr = plots_set_ref()$drugs$sim_vs_control,
          sim_vs_all = plots_set_ref()$drugs$sim_vs_all,
          community = plots_set_ref()$drugs$community,
          control = plots_set_ref()$drugs$control,
          drug = plots_set_ref()$drugs$drug,
          clique_ids = plots_set_ref()$drugs$clique_ids
        ) # end data.frame 
        
        # highlight control nodes
        nodes[which(nodes$control == 1),]$color.border <- "#008744"
        nodes[which(nodes$control == 1),]$borderWidth <- 1.5
        
        # assign node color based on selected metric
        if (input$samples_node_color_metric == "none") {
          
          nodes$color.background <- "#f1f1f1"
          nodes$color.highlight.background <- "#f1f1f1"
          
        } else {
          
          if (input$samples_node_color_metric == "degree")
            col_values <- plots_set_ref()$drugs$degree_samples
          else if (input$samples_node_color_metric == "cells")
            col_values <- plots_set_ref()$drugs$live_cells
          else if (input$samples_node_color_metric == "mfi")
            col_values <- plots_set_ref()$mfis[,input$samples_node_color_channel]
          else if (input$samples_node_color_metric == "sim_vs_ctr")
            col_values <- plots_set_ref()$drugs$sim_vs_control
          else if (input$samples_node_color_metric == "sim_vs_all")
            col_values <- plots_set_ref()$drugs$sim_vs_all
          
          col_fun <- colorRamp2(c(min(col_values), median(col_values), max(col_values)),
                                c("#0057e7", "#fdf498", "#d62d20"))
          nodes$color.background <- col_fun(col_values)
          nodes$color.highlight.background <- col_fun(col_values)
          
        } # end else
        
        # assign node label based on user input
        if (input$samples_node_label == "well") {
          nodes$label <- plots_set_ref()$drugs$file
        } else if (input$samples_node_label == "drug") {
          
          if ("concentration" %in% colnames(plots_set_ref()$drugs))
            nodes$label <- paste(plots_set_ref()$drugs$drug, plots_set_ref()$drugs$concentration, sep = "_")
          else
            nodes$label <- plots_set_ref()$drugs$drug
          
        } else if (input$samples_node_label == "comm") {
          nodes$label <- plots_set_ref()$drugs$community
        } # end esle if
        
        # make edges data frame
        edges <- data.frame(
          from = plots_set_ref()$samples_network$from,
          to = plots_set_ref()$samples_network$to,
          weight = plots_set_ref()$samples_network$weight
        ) # end data.frame
        
        # assign edge colors
        if (input$samples_edge_color_metric == "weight") {
          edges$color <- plots_set_ref()$samples_network$color_weight
        } else {
          edges$color <- "#cbcbcb"
        } # end esle
        
        # make igraph object
        g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
        
        return(g)
        
      }) # end samples_igraph
      
      
      filtered_samples_igraph <- reactive({
        
        g <- samples_igraph()
        
        validate(need(g, message = FALSE))
        
        g <- delete_vertices(g, v = which(V(g)$degree < input$node_degree_range[1]))
        g <- delete_vertices(g, v = which(V(g)$degree > input$node_degree_range[2]))
        
        g <- delete_vertices(g, v = which(V(g)$sim_vs_ctr < input$sim_vs_ctr_range[1]))
        g <- delete_vertices(g, v = which(V(g)$sim_vs_ctr > input$sim_vs_ctr_range[2]))
        
        g <- delete_vertices(g, v = which(V(g)$sim_vs_all < input$sim_vs_all_range[1]))
        g <- delete_vertices(g, v = which(V(g)$sim_vs_all > input$sim_vs_all_range[2]))
        
        clique_choices <- NULL
        
        for (i in 1:length(V(g))) {
          if (!is.na(V(g)$clique_ids[i]))
            clique_choices <- union(clique_choices, strsplit(V(g)$clique_ids[i], split = ",")[[1]])
        } # end for
        
        # update the node id input box
        updateSelectInput(session, "samples_select_by_id", choices = V(g)$name[order(V(g)$name)])
        updateSelectInput(session, "samples_select_by_drug", choices = V(g)$drug[order(V(g)$drug)])
        updateSelectInput(session, "samples_select_by_comm", choices = V(g)$community[order(as.integer(V(g)$community))])
        updateSelectInput(session, "samples_select_by_clique", choices = clique_choices[order(clique_choices)])
        
        return(g)
        
      }) # end filtered_samples_igraph
      
      
      ### render the well clustering network
      output$samples_network <- renderVisNetwork({
        
        g <- filtered_samples_igraph()
        
        validate(need(g, message = FALSE))
        
        # render the network from igraph object
        visIgraph(g, idToLabel = FALSE, layout = input$samples_network_layout, randomSeed = 211) %>%
          visOptions(highlightNearest = input$samples_highlight_nearest, autoResize = TRUE) %>%
          visInteraction(multiselect = TRUE, hideEdgesOnDrag = TRUE)
        
      }) # end output$samples_network
      
      
      ### observer to select nodes by id
      observeEvent(input$samples_select_by_id, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$samples_select_by_id))
          visNetworkProxy(ns("samples_network")) %>% visUnselectAll()
        else
          visNetworkProxy(ns("samples_network")) %>% visSelectNodes(id = input$samples_select_by_id)
        
      }) # end observeEvent
      
      
      ### observer to select nodes by drug
      observeEvent(input$samples_select_by_drug, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$samples_select_by_drug))
          visNetworkProxy(ns("samples_network")) %>% visUnselectAll()
        else {
          g <- samples_igraph()
          visNetworkProxy(ns("samples_network")) %>%
            visSelectNodes(id = V(g)$name[which(V(g)$drug %in% input$samples_select_by_drug)])
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to select nodes by community
      observeEvent(input$samples_select_by_comm, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$samples_select_by_comm))
          visNetworkProxy(ns("samples_network")) %>% visUnselectAll()
        else {
          g <- samples_igraph()
          visNetworkProxy(ns("samples_network")) %>%
            visSelectNodes(id = V(g)$name[which(V(g)$community %in% input$samples_select_by_comm)]) 
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to select nodes by clique
      observeEvent(input$samples_select_by_clique, ignoreNULL = FALSE, {
        
        ns <- session$ns
        
        if (is.null(input$samples_select_by_clique))
          visNetworkProxy(ns("samples_network")) %>% visUnselectAll()
        else {
          g <- samples_igraph()
          
          ids <- NULL
          for (c in input$samples_select_by_clique) {
            ids <- union(ids, strsplit(plots_set_ref()$cliques[which(plots_set_ref()$cliques$ID == c),]$File, split = ",")[[1]])
          } # end for
          
          visNetworkProxy(ns("samples_network")) %>%
            visSelectNodes(id = ids) 
        } # end else
        
      }) # end observeEvent
      
      
      ### observer to add selected nodes form the samples network into the scatter well table
      observe({
        
        input$samples_network_add_selected
        
        ns <- session$ns
        
        visNetworkProxy(ns("samples_network")) %>% visGetSelectedNodes()
        
        validate(need(input$samples_network_selectedNodes, message = FALSE))
        
        selected_wells(list(wells = input$samples_network_selectedNodes))
        
      }) # end observe
      
      return(selected_wells)
      
    } # end function
    
  ) # end moduleServer
  
} # end Vis_samples_network_S


### Server for Scatterplot tab
# @param id String, id of corresponding Vis_scatterplot_UI module
# @param event_plots_data Reactive expression containing a list of input values for scatterplot visualization
Vis_scatterplot_S <- function(id, event_plots_data) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the event_plots_data
      event_plots_data_ref <- reactive({
        validate(need(event_plots_data(), message = FALSE))
        return(event_plots_data())
      }) # end event_plots_data_ref
      
      
      ### return the selected height of the static scatterplot
      scatterplot_height <- function() {
        return(as.integer(input$height))
      } # end scatterplot_height
      
      
      ### function to make the static scatterplot
      output$scatterplot <- renderPlot(height = scatterplot_height, {
        
        d <- event_plots_data_ref()$d
        
        x_lab <- names(event_plots_data_ref()$channel_choices)[which(event_plots_data_ref()$channel_choices == event_plots_data_ref()$scatter_x)]
        y_lab <- names(event_plots_data_ref()$channel_choices)[which(event_plots_data_ref()$channel_choices == event_plots_data_ref()$scatter_y)]
        
        colnames(d)[1] <- "x"
        colnames(d)[2] <- "y"
        
        if (input$color == "none") { # no point colors
          
          p <- ggplot(d, aes(x = x, y = y)) +
            geom_point(size = as.double(input$point_size), alpha = 0.7) +
            labs(x = x_lab, y = y_lab) +
            theme_light() +
            theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20))
          
          if (input$margins == "none")
            p
          else
            ggMarginal(p, type = input$margins, size = 5)
          
        } else if (input$color == "density") { # color by density
          
          p <- ggplot(d, aes(x = x, y = y, color = density)) +
            geom_point(size = as.double(input$point_size), alpha = 0.8) +
            labs(x = x_lab, y = y_lab) +
            theme_light() +
            theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
                  legend.text = element_text(size = 14, angle = 45, hjust = 1), legend.title = element_text(size = 14),
                  axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
            scale_color_viridis_c()
          
          if (input$margins == "none")
            p
          else
            ggMarginal(p, type = input$margins, size = 5)
          
        } else { # color by variables
          
          p <- ggplot(d, aes(x = x, y = y, color = var_name)) +
            geom_point(size = as.double(input$point_size), alpha = 0.8) +
            labs(x = x_lab, y = y_lab) +
            theme_light() +
            theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
                  legend.title = element_blank(), legend.text = element_text(size = 14),
                  axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 5))) +
            scale_color_manual(values = setNames(unique(d$var_color), unique(d$var_name)))
          
          if (input$margins == "none")
            p
          else
            ggMarginal(p, type = input$margins, size = 5, groupColour = TRUE, groupFill = TRUE)
          
        } # end else
        
      }) # end output$scatterplot
      
    } # end function
    
  ) # end moduleServer
  
} # end Vis_scatterplot_S


### Server for Violins tab
# @param id String, id of corresponding Vis_violins_UI module
# @param event_plots_data Reactive expression containing a list of input values for violins visualization
Vis_violins_S <- function(id, event_plots_data) {
  
  moduleServer(
    id,
    
    function(input, output, session) {
      
      ### validate and reference the event_plots_data
      event_plots_data_ref <- reactive({
        validate(need(event_plots_data(), message = FALSE))
        return(event_plots_data())
      }) # end event_plots_data_ref
      
      
      ### update UI elements when new data is received
      updateSelectInput(session, "channels", choices = event_plots_data_ref()$channel_choices, selected = event_plots_data_ref()$channel_choices)
      
      
      ### data for the channel intensities plot
      violins_data <- reactive({
        
        validate(need(input$channels, message = FALSE))
        
        d <- event_plots_data_ref()$data
        
        if (input$color == "none") {
          
          data_long <- NULL
          for (c in input$channels)
            data_long <- rbind(data_long, cbind(d[,c], rep(names(event_plots_data_ref()$channel_choices)[which(event_plots_data_ref()$channel_choices == c)], nrow(d))))
          
          data_long <- as.data.frame(data_long)
          
          colnames(data_long) <- c("intensity", "channel")
          
          data_long$intensity <- as.double(data_long$intensity)
          data_long$channel <- as.factor(data_long$channel)
          
        } else {
          
          data_long <- NULL
          for (c in input$channels)
            data_long <- rbind(data_long, cbind(d[,c], rep(names(event_plots_data_ref()$channel_choices)[which(event_plots_data_ref()$channel_choices == c)], nrow(d)),
                                                d[,event_plots_data_ref()$selected_var],
                                                d[,paste0(event_plots_data_ref()$selected_var, "_color")]))
          
          data_long <- as.data.frame(data_long)
          
          colnames(data_long) <- c("intensity", "channel", "var_name", "var_color")
          
          data_long$intensity <- as.double(data_long$intensity)
          data_long$channel <- as.factor(data_long$channel)
          data_long$var_name <- as.factor(data_long$var_name)
          data_long$var_color <- as.character(data_long$var_color)
          
        } # end else
        
        return(data_long)
        
      }) # end violins_data()
      
      
      ### width of violins plot, depending on selected type, auto (fill width of the browser window) or 2000px
      plot_width <- function() {
        
        if (input$width_type == "auto")
          return("auto")
        else
          return(2000)
        
      } # end plot_width 
      
      
      ### render the channel intensities plot
      output$violins_plot <- renderPlot(width = plot_width, {
        
        validate(need(violins_data(), message = FALSE))
        
        d <- violins_data()
        
        violins_var <- isolate(input$color)
        
        if (violins_var == "none") {
          
          ggplot(d, aes(x = channel, y = intensity)) +
            geom_violin(trim = FALSE, fill = "#f1f1f1", alpha = 0.7) +
            geom_boxplot(width = 0.1, fill = "#f1f1f1", alpha = 0.7) +
            labs(x = NULL) +
            theme_light() +
            theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
          
        } else {
          
          ggplot(d, aes(x = channel, y = intensity, fill = var_name)) +
            geom_violin(trim = FALSE, alpha = 0.7) +
            geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), alpha = 0.7) +
            labs(x = NULL) +
            theme_light() +
            theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
                  legend.text = element_text(size = 14), legend.title = element_blank()) +
            scale_fill_manual(values = setNames(unique(d$var_color), unique(d$var_name)))
          
        } # end else
        
      }) # end output$violins_plot
      
    } # end function
    
  ) # end moduleServer
  
} # end Vis_violins_S


########################################################################################################################
### --- helper functions --- ###########################################################################################
########################################################################################################################

### read COMPARE output and FCS files and generate data to be visualized in Plots tab
# @param channels Character vector of channel names from channels_step9 input
# @param paths Reactive expression containing a list of paths to data and out folders
# @return plots_set List of data frames for visualizations in Plots view
get_plots_set <- function(channels, paths) {
  
  # load the COMPARE RData file from Step 8
  e = new.env()
  rdata_obj <- load(paste0(paths()$out, "compare_clustering.RData"), envir = e)
  rdata_obj <- e[[rdata_obj]]
  
  ### drugs table
  drugs <- read.table(paste0(paths()$out, "drugs_table.tsv"), header = TRUE, sep = "\t")
  
  drugs$control <- as.character(drugs$control)
  drugs$community <- as.character(drugs$community)
  
  
  #### make data frame for samples network
  
  samples_network <- as.data.frame(cbind(as_edgelist(rdata_obj$samples_graph, names = TRUE),
                                             as.double(E(rdata_obj$samples_graph)$weight)))
  colnames(samples_network) <- c("from","to","weight")
  
  samples_network$weight <- as.double(samples_network$weight)
  
  # option to color edges by weight for when the network is rendered
  col_fun <- colorRamp2(c(min(samples_network$weight), median(samples_network$weight), max(samples_network$weight)),
                        c("#008744", "#f1f1f1", "#a200ff"))
  samples_network$color_weight <- col_fun(samples_network$weight)
  
  # round similarity scores
  drugs$sim_vs_control <- round(drugs$sim_vs_control, 3)
  drugs$sim_vs_all <- round(drugs$sim_vs_all, 3)
  
  # calculate node degrees in samples network
  drugs$degree_samples <- NA
  for (i in 1:nrow(drugs)) {
    drugs[i,]$degree_samples <- length(which(samples_network$from == drugs[i,]$file)) + length(which(samples_network$to == drugs[i,]$file))
  } # end for
  
  
  ### make data frame for dispersion network
  
  disp_network <- as.data.frame(cbind(as_edgelist(rdata_obj$dispersion_graph, names = TRUE),
                                          as.double(E(rdata_obj$dispersion_graph)$weight)))
  colnames(disp_network) <- c("from","to","weight")
  
  disp_network$weight <- as.double(disp_network$weight)
  
  # calculate node degrees in samples network
  drugs$degree_disp <- NA
  for (i in 1:nrow(drugs)) {
    drugs[i,]$degree_disp <- length(which(disp_network$from == drugs[i,]$file)) + length(which(disp_network$to == drugs[i,]$file))
  } # end for
  
  # option to color edges by weight for when the network is rendered
  col_fun <- colorRamp2(c(min(disp_network$weight), median(disp_network$weight), max(disp_network$weight)),
                        c("#008744", "#f1f1f1", "#a200ff"))
  disp_network$color_weight <- col_fun(disp_network$weight)
  
  
  ### read the FCS files and get the number of cells, MFI values, and channels table
  
  mfis <- matrix(data = NA, nrow = nrow(drugs), ncol = length(channels) + 1)
  colnames(mfis) <- c("well", channels)
  mfis[,"well"] <- drugs$file
  
  for (i in 1:nrow(drugs)) {
    
    well <- drugs[i,]$file
    fname <- paste0(well, ".fcs")
    
    # read FCS file
    ff <- read.FCS(paste0(paths()$data, fname), transformation = FALSE)
    
    # channels table
    channels_table <- ff@parameters@data[,c("name","desc")]
    
    # get MFI values for provided channel names
    expr_data <- ff@exprs[,channels, drop = FALSE]
    
    for (channel in channels) {
      
      dt_tmp <- expr_data[,channel]
      
      ### Morteza's code to calculate MFI
      dt_tmp = dt_tmp[which(0 <= dt_tmp)]                         # N.B.: which igonores NA and NaNs. Non-positives are always non-positives even in presence of drift
      IQR_ = IQR(dt_tmp)                                          # inter-quantile region
      quartiles_ = quantile(dt_tmp, probs = c(.25, .75))          # 25th and 75th percentiles
      lowWhisker_ = max(min(dt_tmp), quartiles_[1] - IQR_*1.5)    # lower whisker
      upWhisker_ = min(max(dt_tmp), quartiles_[2] + IQR_*1.5)
      mfis[which(mfis[,"well"] == well), channel] = round(median(dt_tmp[lowWhisker_ < dt_tmp & dt_tmp < upWhisker_]), 3)
      
    } # end for
    
  } # end for
  
  mfis <- as.data.frame(mfis)
  channels_table <- as.data.frame(channels_table)
  
  rownames(channels_table) <- channels_table$name
  
  # mark which plots were processed for Plots tabs
  channels_table$processed_for_plots <- NA
  channels_table[channels,]$processed_for_plots <- "yes"
  
  
  colnames(mfis) <- c("well", channels_table[channels,]$desc)
  
  for (desc in channels_table[channels,]$desc) {
    mfis[,desc] <- as.double(mfis[,desc])
  } # end for
  
  
  ####### cliques UMAP and centroids
  
  umap <- rdata_obj$umap_
  centroids <- read.table(paste0(paths()$out, "Centroids.tsv"), sep = "\t", header = TRUE, row.names = 1)
  
  colnames(centroids) <- channels_table[channels,]$desc
  
  
  ####### cliques table
  
  cliques <- read.table(paste0(paths()$out, "Cliques.tsv"), sep = "\t", header = TRUE)
  
  cliques$ID <- paste0("C", cliques$ID) # to match cliques ID in umap and centroids tables
  cliques$length <- NA # number of wells in clique
  cliques$total_live_cells <- NA # total events in clique
  
  for (i in 1:nrow(cliques)) {
    cliques[i,]$length <- length(strsplit(cliques[i,]$File, split = ",")[[1]])
    cliques[i,]$total_live_cells <- sum(as.numeric(strsplit(cliques[i,]$live_cels, split = ",")[[1]]))
    
  } # end for 
  
  drugs$clique_ids <- NA
  
  for (i in 1:nrow(drugs)) {
    
    for (c in 1:nrow(cliques)) {
      
      if (drugs$file[i] %in% strsplit(cliques[c,]$File, split = ",")[[1]]) {
        
        if (is.na(drugs[i,]$clique_ids))
          drugs[i,]$clique_ids <- cliques[c,]$ID
        else
          drugs[i,]$clique_ids <- paste(c(drugs[i,]$clique_ids, cliques[c,]$ID), collapse = ",")
        
      } # end if
      
    } # end for
    
  } # end for
  
  # separate control cluster wells from the rest
  ctr_wells <- drugs[which(drugs$community == "0"),]
  
  
  return(
    list(
      drugs = drugs,
      samples_network = samples_network,
      disp_network = disp_network,
      mfis = mfis,
      channels = channels,
      channels_table = channels_table,
      umap = umap,
      centroids = centroids,
      cliques = cliques,
      ctr_wells = ctr_wells
    ) # end list
  ) # end return
  
} # end get_plots_set


### Function to get density of points in 2 dimensions, source: https://slowkow.com/notes/ggplot2-color-by-density/
# @param x Numeric vector
# @param y Numeric vector
# @param n Create a square n by n grid to compute density
# @return The density within each square
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
} # end get_density


########################################################################################################################
### --- main Sniny app UI --- ##########################################################################################
########################################################################################################################

### generate UI for Workflow and Tables views and call plots_UI() to generate UI for Plots view
ui <- navbarPage(
  
  title = "COMPARE-suite",
  
  theme = shinytheme("simplex"),
  
  tabPanel(
    title = "Workflow",
    useShinyjs(),
    
    # change body background
    tags$head(tags$style(HTML("body {background-color:#ffffff;}"))),
    
    # change background of active tabpanel title
    tags$style(HTML(".tabbable > .nav > li[class=active] > a {background-color:#ffffff;}")),
    
    fluidRow(
      
      column(
        8,
        
        fluidRow(
          
          column(
            6,
            div(textInput("data_path", "Data path:", placeholder = "Path to Annotations.txt and FCS files", value = paste0("..", .Platform$file.sep, "data", .Platform$file.sep), width = "100%"),
                style = "font-size:95%"),
            bsPopover("data_path", title = NULL,
                      content = "Path to folder containing Annotations.txt and FCS files.</br><b>Warning:</b> Step 1 rewrites Annotations.txt file and Steps 2 to 4 rewrite the FCS files. We recommend to have backup copies of the annotation and FCS files.",
                      placement = "bottom", trigger = "focus", options = list(container = "body"))
          ), # end column
          
          column(
            6,
            div(textInput("out_path", "Out path:", placeholder = "Path to output folder", value = paste0("..", .Platform$file.sep, "out", .Platform$file.sep), width = "100%"),
                style = "font-size:95%"),
            bsPopover("out_path", title = NULL, content = "Path to folder for output files",
                      placement = "bottom", trigger = "focus", options = list(container = "body"))
          ) # end column
          
        ), # end fluidRow
        
        ### STEP 1
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step1", "1 - Import", value = FALSE),
                style = "font-size:113%; color:#333333"),
          ), # end column
          column(
            8,
            fluidRow(
              column(
                4,
                div(uiOutput("uio_min_events"), style = "font-size:95%"),
                bsPopover("uio_min_events", title = NULL, content = "Wells with events below this number are skipped",
                          placement = "right", trigger = "focus", options = list(container = "body"))
              ) # end column
            ) # end fluidRow
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### STEP 2
        div(checkboxInput("include_step2", "2 - Spillover compensation", value = FALSE),
            style = "font-size:113%; color:#333333"),
        
        br(),
        
        
        ### STEP 3
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step3", "3 - Signal drift correction", value = FALSE),
                style = "font-size:113%; color:#333333")
          ), # end column
          column(
            8,
            div(uiOutput("uio_channels_step3"), style = "font-size:95%"),
            bsPopover("uio_channels_step3", title = NULL, content = "Channels to use for signal drift correction",
                      placement = "right", trigger = "focus", options = list(container = "body")),
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_drctn_step3"), style = "font-size:95%"),
                bsPopover("uio_drctn_step3", title = NULL, content = "Direction of signal drift",
                          placement = "right", trigger = "focus", options = list(container = "body")),
              ) # end column
            ), # end fluidRow
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_correct_step3"), style = "font-size:95%")
              ), # end column
              column(
                4,
                div(uiOutput("uio_fitplot_step3"), style = "font-size:95%")
              ), # end column
              column(
                4,
                div(uiOutput("uio_heatplot_step3"), style = "font-size:95%")
              ) # end column
            ) # end fluidRow
            
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### STEP 4
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step4", "4 - Viability correction", value = FALSE),
                style = "font-size:113%; color:#333333")
          ), # end column
          column(
            8,
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_drctn_step4"), style = "font-size:95%"),
                bsPopover("uio_drctn_step4", title = NULL, content = "Direction of viability bias",
                          placement = "right", trigger = "focus", options = list(container = "body")),
              ) # end column
            ), # end fluidRow
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_correct_step4"), style = "font-size:95%")
              ), # end column
              column(
                4,
                div(uiOutput("uio_fitplot_step4"), style = "font-size:95%")
              ), # end column
              column(
                4,
                div(uiOutput("uio_heatplot_step4"), style = "font-size:95%")
              ) # end column
            ) # end fluidRow
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### STEP 5
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step5", "5 - Similarity matrix generator", value = FALSE),
                style = "font-size:113%; color:#333333")
          ), # end column
          column(
            8,
            div(uiOutput("uio_channels_step5"), style = "font-size:95%"),
            bsPopover("uio_channels_step5", title = NULL, content = "Channels to use for similarity matrix generation",
                      placement = "right", trigger = "focus", options = list(container = "body")),
            fluidRow(
              column(
                4,
                div(uiOutput("uio_n"), style = "font-size:95%"),
                bsPopover("uio_n", title = NULL, content = "Number of subspaces to devide the original space to",
                          placement = "right", trigger = "focus", options = list(container = "body"))
              ), # end column
              column(
                4,
                div(uiOutput("uio_num_cores"), style = "font-size:95%"),
                bsPopover("uio_num_cores", title = NULL, content = "Number of CPU cores to use",
                          placement = "right", trigger = "focus", options = list(container = "body"))
              ) # end column
            ) # end fluidRow
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### STEP 6
        div(checkboxInput("include_step6", "6 - Similarity matrix heatmap", value = FALSE),
            style = "font-size:113%; color:#333333"),
        
        br(),
        
        
        ### STEP 7
        div(checkboxInput("include_step7", "7 - Negative control outlier detector", value = FALSE),
            style = "font-size:113%; color:#333333"),
        
        br(),
        
        
        ### STEP 8
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step8", "8 - Clustering", value = FALSE),
                style = "font-size:113%; color:#333333")
          ), # end column
          column(
            8,
            div(uiOutput("uio_channels_step8"), style = "font-size:95%"),
            bsPopover("uio_channels_step8", title = NULL, content = "Channels to use for clustering",
                      placement = "right", trigger = "focus", options = list(container = "body")),
            fluidRow(
              column(
                4,
                div(uiOutput("uio_nn"), style = "font-size:95%"),
                bsPopover("uio_nn", title = NULL, content = "Number of nearest neighbors for UMAP calculation",
                          placement = "right", trigger = "focus", options = list(container = "body"))
              ) # end column
            ) # end fluidRow
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### STEP 9
        fluidRow(
          column(
            4,
            div(checkboxInput("include_step9", "9 - Make interactive plots", value = FALSE),
                style = "font-size:113%; color:#333333")
          ), # end column
          column(
            8,
            div(uiOutput("uio_channels_step9"), style = "font-size:95%"),
            bsPopover("uio_channels_step9", title = NULL, content = "Channels to include in interactive plots",
                      placement = "right", trigger = "focus", options = list(container = "body"))
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### RUN WORKFLOW
        fluidRow(
          column(
            3,
            actionButton("run_workflow", "Run", width = "100%")
          ), # end column
          column(
            9,
            div(textOutput("run_error_message"), style = "margin-top:9px; color:#c73824")
          ) # end column
        ), # end fluidRow
        
        textOutput("message_output"),
        
        tags$head(
          
          tags$style(
            "#message_output{font-size:12px; background: #ffffff;
                             border: 1px; border-style: solid; border-color: #dddddd; border-radius: 4px;
                             padding: 10px; height: 600px; overflow: auto; margin-top: 20px;}"
          ), # end tags$style
          
          # Shiny to javascript binding to scroll message_output to bottom once called
          # source: https://stackoverflow.com/questions/36677726/r-how-do-i-automatically-scroll-to-the-bottom-of-a-div-in-shinyapp
          tags$script(
            'Shiny.addCustomMessageHandler("scrollCallback",
                  function(x) {
                    var objDiv = document.getElementById("message_output");
                    objDiv.scrollTop = objDiv.scrollHeight;});'
          ) # end tags$script
          
        ) # end tags$head
        
      ) # end column
      
    ) # end fluidRow
    
  ), # end tabPanel
  
  
  tabPanel(
    title = "Plots",
    useShinyjs(),
    
    plots_UI("Plots_UI_1")
    
  ), # end tabPanel
  
  
  tabPanel(
    title = "Tables",
    
    tabsetPanel(
      
      tabPanel(
        title = "Drugs",
        br(),
        div(dataTableOutput("drugs_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Channels",
        br(),
        div(dataTableOutput("channels_table_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Channel MFIs",
        br(),
        div(dataTableOutput("mfis_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques",
        br(),
        div(dataTableOutput("cliques_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Samples edges",
        br(),
        div(dataTableOutput("samples_network_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Dispersion edges",
        br(),
        div(dataTableOutput("disp_network_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques UMAP",
        br(),
        div(dataTableOutput("umap_dt"), style = "font-size:95%; overflow: auto")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques centroids",
        br(),
        div(dataTableOutput("centroids_dt"), style = "font-size:95%; overflow: auto")
      ) # end tabPanel
      
    ) # end tabsetPanel
    
  ), # end tabPanel
  
  br(),
  br()
) # end ui


########################################################################################################################
### --- main Shiny app Server --- ######################################################################################
########################################################################################################################

### perform server processing for Workflow and Tables views and call plots_S() for server processing for Plots view
server <- function(input, output, session) { 
  
  ### set options
  options(shiny.maxRequestSize = 3000*1024^2,
          stringsAsFactors = FALSE
  ) # end options
  
  ############################################################################################################
  ### --- render UI elements for selected steps --- ##########################################################
  ############################################################################################################
  
  ### render UI elements for additional parameters in step 1
  observeEvent(input$include_step1, {
    
    if (input$include_step1) {
      
      output$uio_min_events <- renderUI({
        numericInput("min_events", "Minimum events:", value = 1000, min = 100, max = 100000, step = 1000, width = "100%")
      }) # end output$uio_min_events
      
    } else {
      
      output$uio_min_events <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent  
    
  
  ### render UI elements for additional parameters in step 3
  observeEvent(input$include_step3, {
    
    if (input$include_step3) {
      
      output$uio_channels_step3 <- renderUI({
        div(textInput("channels_step3", "Channels:", width = "100%",
                      placeholder = "Coma-separated names, e.g., SSC-H,VL1-H,VL6-H"),
            style = "font-size:95%")
      }) # end out$uio_channels_step3
      
      output$uio_correct_step3 <- renderUI({
        checkboxInput("correct_step3", "Perform correction", value = TRUE)
      }) # end output$uio_correct_step3
      
      output$uio_drctn_step3 <- renderUI({
        div(selectInput("drctn_step3", "Direction of bias:", choices = c("Column" = "column", "Row" = "row"),
                        width = "100%"),
            style = "font-size:95%")
      }) # end out$uio_channels_step3
      
      output$uio_fitplot_step3 <- renderUI({
        checkboxInput("fitplot_step3", "Plot regressed line", value = TRUE)
      }) # end output$uio_fitplot_step3
      
      output$uio_heatplot_step3 <- renderUI({
        checkboxInput("heatplot_step3", "Plot plate heatmaps", value = TRUE)
      }) # end output$uio_heatplot_step3
      
    } else {
      
      output$uio_channels_step3 <- renderUI({NULL})
      output$uio_correct_step3 <- renderUI({NULL})
      output$uio_drctn_step3 <- renderUI({NULL})
      output$uio_fitplot_step3 <- renderUI({NULL})
      output$uio_heatplot_step3 <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent  
  
  
  ### render UI elements for additional parameters in step 4
  observeEvent(input$include_step4, {
    
    if (input$include_step4) {
      
      output$uio_correct_step4 <- renderUI({
        checkboxInput("correct_step4", "Perform correction", value = TRUE)
      }) # end output$uio_correct_step4
      
      output$uio_drctn_step4 <- renderUI({
        div(selectInput("drctn_step4", "Direction of bias:", choices = c("Column" = "column", "Row" = "row"),
                        width = "100%"),
            style = "font-size:95%")
      }) # end out$uio_channels_step4
      
      output$uio_fitplot_step4 <- renderUI({
        checkboxInput("fitplot_step4", "Plot regressed line", value = TRUE)
      }) # end output$uio_fitplot_step3
      
      output$uio_heatplot_step4 <- renderUI({
        checkboxInput("heatplot_step4", "Plot plate heatmaps", value = TRUE)
      }) # end output$uio_heatplot_step3
      
    } else {
      
      output$uio_correct_step4 <- renderUI({NULL})
      output$uio_drctn_step4 <- renderUI({NULL})
      output$uio_fitplot_step4 <- renderUI({NULL})
      output$uio_heatplot_step4 <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent
  
  
  ### render UI elements for additional parameters in step 5
  observeEvent(input$include_step5, {
    
    if (input$include_step5) {
      
      output$uio_channels_step5 <- renderUI({
        div(textInput("channels_step5", "Channels:", width = "100%",
                      placeholder = "Coma-separated names, e.g., SSC-H,VL1-H,VL6-H"),
            style = "font-size:95%")
      }) # end out$uio_channels_step5
      
      output$uio_n <- renderUI({
        numericInput("n", "Number of subspaces:", value = 3, min = 10, max = 10, width = "100%")
      }) # end output$uio_n
      
      output$uio_num_cores <- renderUI({
        numericInput("num_cores", "Number of cores:", value = 1, min = 1, max = 100, width = "100%")
      }) # end output$uio_num_cores
      
    } else {
      
      output$uio_channels_step5 <- renderUI({NULL})
      output$uio_n <- renderUI({NULL})
      output$uio_num_cores <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent
  
  
  ### render UI elements for additional parameters in step 8
  observeEvent(input$include_step8, {
    
    if (input$include_step8) {
      
      output$uio_channels_step8 <- renderUI({
        div(textInput("channels_step8", "Channels:", width = "100%",
                      placeholder = "Coma-separated names, e.g., SSC-H,VL1-H,VL6-H"),
            style = "font-size:95%")
      }) # end out$uio_channels_step8
      
      output$uio_nn <- renderUI({
        numericInput("nn", "Nearest neighbors:", value = 5, min = 2, max = 100, width = "100%")
      }) # end output$uio_n
      
    } else {
      
      output$uio_channels_step8 <- renderUI({NULL})
      output$uio_nn <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent
  
  
  ### render UI elements for additional parameters in step 9
  observeEvent(input$include_step9, {
    
    if (input$include_step9) {
      
      output$uio_channels_step9 <- renderUI({
        div(textInput("channels_step9", "Channels:", width = "100%",
                      placeholder = "Coma-separated names, e.g., SSC-H,VL1-H,VL6-H"),
            style = "font-size:95%")
      }) # end out$uio_channels_step9
      
    } else {
      
      output$uio_channels_step9 <- renderUI({NULL})
      
    } # end else
    
  }) # end observeEvent
  
  
  ############################################################################################################
  ### --- run COMPARE workflow --- ###########################################################################
  ############################################################################################################
  
  # reactive value to hold the list of data tables for visualizations
  plots_set <- reactiveVal(NULL)
  
  ### validate and reference the data and output paths
  paths <- reactive({
    
    validate(need(input$out_path, message = FALSE))
    validate(need(input$data_path, message = FALSE))
    
    if (!(endsWith(input$out_path, .Platform$file.sep)))
      out_path <- paste0(input$out_path, .Platform$file.sep)
    else
      out_path <- input$out_path
    
    if (!(endsWith(input$data_path, .Platform$file.sep)))
      data_path <- paste0(input$data_path, .Platform$file.sep)
    else
      data_path <- input$data_path
    
    return(list(data = data_path, out = out_path))
    
  }) # end paths
  
  
  ### run workflow steps
  observeEvent(input$run_workflow, {
    
    # check if at least one file is selected
    if (!(TRUE %in% c(input$include_step1,
                      input$include_step2,
                      input$include_step3,
                      input$include_step4,
                      input$include_step5,
                      input$include_step6,
                      input$include_step7,
                      input$include_step8,
                      input$include_step9))) {
      html(id = "run_error_message",
           html = "Please check at least one workflow step to run",
           add = FALSE)
      return(NULL)
    } # end if
    
    # check if all required columns are present in the annotation file
    annot <- read.table(paste0(paths()$data, .Platform$file.sep, "Annotations.txt"), header = TRUE, sep = "\t")
    if (!all(c("file","drug","plate","row","column","concentration","control") %in% colnames(annot))) {
      
      html(id = "run_error_message",
           html = "Error reading annotation, please see below",
           add = FALSE)
      
      html(id = "message_output",
           html = "Please make sure that the annotation file contains the required columns and the column names are:</br>
                   <strong>file</strong>: name of FCS file without extension</br>
                   <strong>drug</strong>: compound name</br>
                   <strong>plate</strong>: plate number</br>
                   <strong>row</strong>: row coordinate, e.g, A, B, C</br>
                   <strong>column</strong>: column coordinate, e.g, 1, 2, 3</br>
                   <strong>concentration</strong>: compound concentration</br>
                   <strong>control</strong>: either 0 or 1 where 0 = drug compound and 1 = negative control compound
                   </br></br>",
           add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      return(NULL)
      
    } # end else
    
    html(id = "run_error_message",
         html = "",
         add = FALSE)
    
    
    ### RUN STEP 1
    if (input$include_step1) {
      
      validate(need(input$min_events, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 1 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step1_import(min_events = input$min_events,
                    inURL = paths()$data)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 1 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 2
    if (input$include_step2) {
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 2 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step2_spillover_compensation(inURL = paths()$data)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 2 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 3
    if (input$include_step3) {
      
      validate(need(strsplit(input$channels_step3, split = '[,]')[[1]], message = FALSE))
      validate(need(input$correct_step3, message = FALSE))
      validate(need(input$drctn_step3, message = FALSE))
      validate(need(input$fitplot_step3, message = FALSE))
      validate(need(input$heatplot_step3, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 3 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step3_signal_drift_correction(chnls_ = strsplit(input$channels_step3, split = '[,]')[[1]],
                                      CORRECT = input$correct_step3,
                                      drctn_ = input$drctn_step3,
                                      FITPLOT = input$fitplot_step3,
                                      HEATPLOT = input$heatplot_step3,
                                      inURL = paths()$data,
                                      outURL = paths()$out)
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 3 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 4
    if (input$include_step4) {
      
      validate(need(input$correct_step4, message = FALSE))
      validate(need(input$drctn_step4, message = FALSE))
      validate(need(input$fitplot_step4, message = FALSE))
      validate(need(input$heatplot_step4, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 4 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step4_viability_correction(CORRECT = input$correct_step4,
                                   drctn_ = input$drctn_step4,
                                   FITPLOT = input$fitplot_step4,
                                   HEATPLOT = input$heatplot_step4,
                                   inURL = paths()$data,
                                   outURL = paths()$out)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 4 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 5
    if (input$include_step5) {
      
      validate(need(strsplit(input$channels_step5, split = '[,]')[[1]], message = FALSE))
      validate(need(input$n, message = FALSE))
      validate(need(input$num_cores, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 5...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step5_similarity_matrix_generator(chnls_ = strsplit(input$channels_step5, split = '[,]')[[1]],
                                          n_ = input$n,
                                          cor_ = input$num_cores,
                                          inURL = paths()$data,
                                          outURL = paths()$out)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 5 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 6
    if (input$include_step6) {
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 6 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step6_similarity_matrix_heatmap(outURL = input$out_path)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 6 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 7
    if (input$include_step7) {
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 7 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step7_negative_control_outlier_detector(inURL = paths()$data,
                                                outURL = paths()$out)
        
      },
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 7 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 8
    if (input$include_step8) {
      
      validate(need(strsplit(input$channels_step8, split = '[,]')[[1]], message = FALSE))
      validate(need(input$nn, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 8 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # run with call handlers to redirect messages
      withCallingHandlers({
        
        step8_clustering(chnls_ = strsplit(input$channels_step8, split = '[,]')[[1]],
                         nn_ = input$nn,
                         inURL = paths()$data,
                         outURL = paths()$out)
        
      }, 
      
      message = function(m) {
        
        # print message
        html(id = "message_output", html = paste0(Sys.time(), " : ", m$message, "</br>"), add = TRUE)
        
        # scroll window
        session$sendCustomMessage(type = "scrollCallback", 1)
        
      } # end message
      
      ) # end withCallingHandlers
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 8 done!</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
    
    ### RUN STEP 9
    if (input$include_step9) {
      
      validate(need(input$channels_step9, message = FALSE))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Running Step 9 ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # list of files in output directory
      files <- list.files(paths()$out)
      
      # check if the needed files are present
      
      if (!("compare_clustering.RData" %in% files)) {
        html(id = "run_error_message",
             html = "Can't find compare_clustering.RData file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      if (!("drugs_table.tsv" %in% files)) {
        html(id = "run_error_message",
             html = "Can't find drugs_table.tsv file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      if (!("Cliques.tsv" %in% files)) {
        html(id = "run_error_message",
             html = "Can't find Cliques.tsv file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      html(id = "message_output", html = paste0(Sys.time(), " : Reading COMPARE output and FCS files, please wait ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      channels <- strsplit(input$channels_step9, split = '[,]')[[1]]
      
      ### generate data for interactive plots
      plots_set(get_plots_set(channels, paths))
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 9 done! Proceed to Plots tab</br>"), add = TRUE)
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
  }) # end observeEvent
  
  
  ######################################################################################################################
  ### --- Plots view --- ###############################################################################################
  ######################################################################################################################
  
  plots_S("Plots_UI_1", plots_set, paths)
  
  
  
  ######################################################################################################################
  ### --- Tables view --- ##############################################################################################
  ######################################################################################################################
  
  ### render the drugs table
  output$drugs_dt <- renderDT(
    plots_set()$drugs[,c("file","drug","concentration","live_cells","control","community","sim_vs_control","sim_vs_all","clique_ids")],
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$drugs_dt
  
  
  ### render the channels table
  output$channels_table_dt <- renderDT(
    plots_set()$channels_table,
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$channels_table_dt
  
  ### render the channel MFIs table
  output$mfis_dt <- renderDT(
    plots_set()$mfis,
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$mfis_dt
  
  ### render the cliques table
  output$cliques_dt <- renderDT(
    plots_set()$cliques,
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$cliques_dt
  
  ### render the samples network table
  output$samples_network_dt <- renderDT(
    plots_set()$samples_network[,c("from","to","weight")],
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$samples_network_dt
  
  ### render the samples network table
  output$disp_network_dt <- renderDT(
    plots_set()$disp_network[,c("from","to","weight")],
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$disp_network_dt
  
  ### render the cliques UMAP table
  output$umap_dt <- renderDT(
    plots_set()$umap,
    selection = "none",
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$umap_dt
  
  ### render the cliques centroids table
  output$centroids_dt <- renderDT(
    plots_set()$centroids,
    selection = "none",
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$umap_dt
  
  
} # end server

### make Shiny app
shinyApp(ui = ui, server = server)
