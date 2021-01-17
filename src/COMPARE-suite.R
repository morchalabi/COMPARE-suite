
require(shiny)
require(shinythemes)
require(shinyjs)
require(shinyFiles)
require(shinyBS)
require(plotly)
require(ggExtra)
require(DT)
require(visNetwork)
require(igraph)
require(circlize)
require(inlmisc)
require(ComplexHeatmap)
require(writexl)

# load the workflow functions into global environment
source("functions/1_import.R")
source("functions/2_spillover_compensation.R")
source("functions/3_signal_drift_correction.R")
source("functions/4_viability_correction.R")
source("functions/5_similarity_matrix_generator.R")
source("functions/6_similarity_matrix_heatmap.R")
source("functions/7_negative_control_outlier_detector.R")
source("functions/8_clustering.R")

# Function to get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
#
# source: https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
} # end get_density

########################################################################################################################
########################################################################################################################
### --- Shiny UI --- ###################################################################################################
########################################################################################################################
########################################################################################################################

ui <- navbarPage(
  
  title = "COMPARE-suite",
  
  theme = shinytheme("simplex"),
  
  tabPanel(
    title = "Workflow", ################################################################################################
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
                      content = "Relative path to folder containing Annotations.txt and FCS files.</br><b>Warning:</b> Step 1 rewrites Annotations.txt file and Steps 2 to 4 rewrite the FCS files. We recommend to have backup copies of the annotation and FCS files.",
                      placement = "bottom", trigger = "hover")
          ), # end column
          
          column(
            6,
            div(textInput("out_path", "Out path:", placeholder = "Path to output folder", value = paste0("..", .Platform$file.sep, "out", .Platform$file.sep), width = "100%"),
                style = "font-size:95%"),
            bsPopover("out_path", title = NULL, content = "Relative path to folder for output files",
                      placement = "bottom", trigger = "hover")
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
                          placement = "bottom", trigger = "hover")
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
                      placement = "bottom", trigger = "hover"),
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_drctn_step3"), style = "font-size:95%"),
                bsPopover("uio_drctn_step3", title = NULL, content = "Direction of signal drift",
                          placement = "right", trigger = "hover"),
              ) # end column
            ), # end fluidRow
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_correct_step3"), style = "font-size:95%")
                #bsPopover("uio_correct_step3", title = NULL, content = "Whether to perform signal drift correction",
                #          placement = "bottom", trigger = "hover")
              ), # end column
              column(
                4,
                div(uiOutput("uio_fitplot_step3"), style = "font-size:95%")
                #bsPopover("uio_fitplot_step3", title = NULL, content = "Whether to plot regressed lines",
                #          placement = "bottom", trigger = "hover")
              ), # end column
              column(
                4,
                div(uiOutput("uio_heatplot_step3"), style = "font-size:95%")
                #bsPopover("uio_heatplot_step3", title = NULL, content = "Whether to plot plate heatmaps",
                #          placement = "bottom", trigger = "hover")
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
                          placement = "right", trigger = "hover"),
              ) # end column
            ), # end fluidRow
            
            fluidRow(
              column(
                4,
                div(uiOutput("uio_correct_step4"), style = "font-size:95%")
                #bsPopover("uio_correct_step4", title = NULL, content = "Whether to perform cell viability correction",
                #          placement = "bottom", trigger = "hover")
              ), # end column
              column(
                4,
                div(uiOutput("uio_fitplot_step4"), style = "font-size:95%")
                #bsPopover("uio_fitplot_step4", title = NULL, content = "Whether to plot regressed lines",
                #          placement = "bottom", trigger = "hover")
              ), # end column
              column(
                4,
                div(uiOutput("uio_heatplot_step4"), style = "font-size:95%")
                #bsPopover("uio_heatplot_step4", title = NULL, content = "Whether to plot plate heatmaps",
                #          placement = "bottom", trigger = "hover")
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
                      placement = "bottom", trigger = "hover"),
            fluidRow(
              column(
                4,
                div(uiOutput("uio_n"), style = "font-size:95%"),
                bsPopover("uio_n", title = NULL, content = "Number of subspaces to devide the original space to",
                          placement = "bottom", trigger = "hover")
              ), # end column
              column(
                4,
                div(uiOutput("uio_num_cores"), style = "font-size:95%"),
                bsPopover("uio_num_cores", title = NULL, content = "Number of CPU cores to use",
                          placement = "bottom", trigger = "hover")
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
                      placement = "bottom", trigger = "hover"),
            fluidRow(
              column(
                4,
                div(uiOutput("uio_nn"), style = "font-size:95%"),
                bsPopover("uio_nn", title = NULL, content = "Number of nearest neighbors for UMAP calculation",
                          placement = "bottom", trigger = "hover")
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
            bsPopover("uio_channels_step9", title = NULL, content = "Channels to include in interactive visualizations",
                      placement = "bottom", trigger = "hover")
          ) # end column
        ), # end fluidRow
        
        br(),
        
        
        ### RUN WORKFLOW
        fluidRow(
          column(
            2,
            actionButton("run_workflow", "Run", width = "100%")
          ), # end column
          column(
            10,
            div(textOutput("run_message_center"), style = "margin-top:9px; color:#c73824")
          ) # end column
        ), # end fluidRow
        
        textOutput("message_output"),
        
        tags$head(
          
          tags$style(
            "#message_output{font-size:12px; background: #ffffff;
                             border: 1px; border-style: solid; border-color: #dddddd; border-radius: 4px;
                             padding: 10px; height: 638px; overflow: auto; margin-top: 22px;}"
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
    title = "Plots", ###################################################################################################
    useShinyjs(),
    
    # change body background
    tags$head(tags$style(HTML("body {background-color:#ffffff;}"))),
    
    # change background of active tabpanel title
    tags$style(HTML(".tabbable > .nav > li[class=active] > a {background-color:#ffffff;}")),
    
    fluidRow(
      
      column(
        8,
        
        tabsetPanel(
          
          tabPanel(
            title = "UMAP",
            br(),
            
            fluidRow(
              
              column(
                4,
                div(selectInput("umap_channel", label = "Channel:", choices = NULL, selected = NULL, width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                2,
                div(selectInput("umap_x", "X axis:", choices = NULL, selected = NULL, width = "100%"), style = "font-size:95%")
              ),
              
              column(
                2,
                div(selectInput("umap_y", "Y axis:", choices = NULL, selected = NULL, width = "100%"), style = "font-size:95%")
              ), # end column
              
              column(
                2,
                div(selectInput("umap_dot_size", "Clique size:", choices = c("None" = "none", "Events" = "total_live_cells", "Wells" = "length"), selected = "total_live_cells", width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                2,
                div(selectInput("umap_labels", "Clique label:", choices = c("None", "Clique ID", "Community"), selected = "Clique ID"),
                    style = "font-size:95%")
              ) # end column
              
            ), # end fluidRow
            
            fluidRow(
              column(
                2,
                textOutput("umap_plot_hover_name"),
                textOutput("umap_plot_hover_centr")
              ), # end column
              column(
                10,
                textOutput("umap_plot_hover_wells"),
                textOutput("umap_plot_hover_drugs")
              ) # end column
            ), # end fluidRow
            
            plotOutput( "umap_plot", height = "760px", click = "umap_plot_click", hover = hoverOpts(id = "umap_plot_hover", delay = 15))
            
          ), # end tabPanel
          
          
          tabPanel(
            title = "Dispersion network",
            br(),
            
            tabsetPanel(
              
              tabPanel(
                title = "Network settings",
                br(),
                
                fluidRow(
                  
                  column(
                    2,
                    div(selectInput("disp_network_layout", "Layout:",
                                    choices = c("FR" = "layout_with_fr", "Hierarchical" = "hierarchical"),
                                    selected = "layout_with_fr", width = "100%"),
                        style = "font-size:95%"),
                    bsPopover("disp_network_layout", title = NULL, content = "FR = Fruchterman-Reingold",
                              placement = "right", trigger = "hover")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("disp_node_color_metric", "Node color:",
                                    choices = c("None" = "none", "MFI" = "mfi", "Degree" = "degree", "Cells" = "cells",
                                                "Sim vs ctr" = "sim_vs_ctr", "Sim vs all" = "sim_vs_all"),
                                    selected = "none"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(selectInput("disp_node_color_channel", label = "Channel for node MFI:",
                                    choices = NULL, selected = NULL, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("disp_node_label", "Node label:",
                                    choices = c("Well" = "well", "Drug" = "drug", "Community" = "comm"),
                                    selected = "well"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("disp_edge_color_metric", "Edge color:",
                                    choices = c("None" = "none", "Weight" = "weight"),
                                    selected = "none"),
                        style = "font-size:95%")
                  ) # end column
                  
                ) # end fluidRow
                
              ), # end tabPanel
              
              tabPanel(
                title = "Select nodes",
                br(),
                
                fluidRow(
                  
                  column(
                    4,
                    div(selectInput("disp_select_by_id", label = "Well ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(selectInput("disp_select_by_drug", label = "Drug:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(selectInput("disp_select_by_comm", label = "Community:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ) # end column
                  
                ) # end fluidRow
                
              ) # end tabPanel
              
            ), # end tabsetPanel
            
            fluidRow(
              
              column(
                10,
                div(checkboxInput("disp_highlight_nearest", "Highlight nearest neighbors", value = FALSE, width = "100%"),
                    style = "font-size:95%; margin-top:-10px")
              ), # end column
              
              column(
                2,
                div(actionButton("disp_network_add_selected", "Add to table" , width = "100%"), style = "margin-bottom: 20px")
              ) # end column
              
            ), # end fluidRow
            
            div(
              style = "background: #ffffff; border: 1px; border-style: solid; border-color: #cccccc; border-radius: 4px;
                       margin-bottom: 22px; padding: 0px",
              visNetworkOutput("disp_network", height = 660)
            ) # end div
            
          ), # end tabPanel
          
          tabPanel(
            title = "Samples network",
            br(),
            
            tabsetPanel(
              
              tabPanel(
                title = "Network settings",
                br(),
                
                fluidRow(
                  
                  column(
                    2,
                    div(selectInput("samples_network_layout", "Layout:",
                                    choices = c("FR" = "layout_with_fr", "KK" = "layout_with_kk"),
                                    selected = "layout_with_fr"),
                        style = "font-size:95%"),
                    bsPopover("samples_network_layout", title = NULL, content = "FR = Fruchterman-Reingold; KK = Kamada-Kawai",
                              placement = "right", trigger = "hover")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("samples_node_color_metric", "Node color:",
                                     choices = c("None" = "none", "MFI" = "mfi", "Degree" = "degree", "Cells" = "cells", "Sim vs ctr" = "sim_vs_ctr", "Sim vs all" = "sim_vs_all"),
                                     selected = "none"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(selectInput("samples_node_color_channel", "Channel for node MFI:", choices = NULL, multiple = FALSE, selected = NULL, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("samples_node_label", "Node label:", choices = c("Well" = "well", "Drug" = "drug", "Community" = "comm"), selected = "well"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    2,
                    div(selectInput("samples_edge_color_metric", "Edge color:", choices = c("None" = "none", "Weight" = "weight"), selected = "none"),
                        style = "font-size:95%")
                  ) # end column
                  
                ) # end fluidRow
                
              ), # end tabPanel
              
              tabPanel(
                title = "Select nodes",
                br(),
                
                fluidRow(
                  
                  column(
                    3,
                    div(selectInput("samples_select_by_id", label = "Well ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    3,
                    div(selectInput("samples_select_by_drug", label = "Drug:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    3,
                    div(selectInput("samples_select_by_comm", label = "Community:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    3,
                    div(selectInput("samples_select_by_clique", label = "Clique ID:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
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
                    div(sliderInput("node_degree_range", "Node degree:", min = 1, max = 100, value = c(1,100), step = 1, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(sliderInput("sim_vs_ctr_range", "Similarity vs control:", min = 0, max = 100, value = c(0,100), step = 0.1, width = "100%"),
                        style = "font-size:95%")
                  ), # end column
                  
                  column(
                    4,
                    div(sliderInput("sim_vs_all_range", "Similarity vs all:", min = 0, max = 100, value = c(0,100), step = 0.1, width = "100%"),
                        style = "font-size:95%")
                  ) # end column
                  
                ) # end fluidRow
                
              ) # end tabPanel
              
            ), # end tabsetPanel
            
            fluidRow(
              
              column(
                10,
                div(checkboxInput("samples_highlight_nearest", "Highlight nearest neighbors", value = FALSE, width = "100%"),
                    style = "font-size:95%; margin-top:-10px")
              ), # end column
              
              column(
                2,
                div(actionButton("samples_network_add_selected", "Add to table" , width = "100%"), style = "margin-bottom: 20px")
              ) # end column
              
            ), # end fluidRow
            
            div(
              style = "background: #ffffff; border: 1px; border-style: solid; border-color: #cccccc; border-radius: 4px;
                       margin-bottom: 22px; padding: 0px",
              visNetworkOutput("samples_network", height = 660)
            ) # end div
            
          ) # end tabPanel
          
        ) # end tabsetPanel
        
      ), # end column
      
      column(
        4,
        
        tabsetPanel(
          
          tabPanel(
            title = "Selected wells",
            br(),
            
            fluidRow(
              column(
                4,
                div(actionButton("scatter_wells_select_all_rows", "Select all", width = "100%"), style = "margin-bottom:20px")
              ), # end column
              column(
                4,
                div(actionButton("scatter_wells_clear_selected_rows", "Clear selected", width = "100%"), style = "margin-bottom:20px")
              ), # end column
              column(
                4,
                div(actionButton("scatter_wells_clear_table", "Clear table", width = "100%"), style = "margin-bottom:20px")
              ) # end column
            ), # end fluidRow
            
            div(
              style = "margin-bottom: 22px; font-size: 89%; height: 720px; overflow: auto;",
              dataTableOutput("scatter_wells_dt")
            ) # end div
            
          ), # end tabPanel
          
          tabPanel(
            title = "Control cluster wells",
            br(),
            
            fluidRow(
              column(
                4,
                div(actionButton("ctr_wells_select_all_rows", "Select all", width = "100%"), style = "margin-bottom:20px")
              ), # end column
              column(
                4,
                div(actionButton("ctr_wells_clear_selected_rows", "Clear selected", width = "100%"), style = "margin-bottom:20px")
              ) # end column
            ), # end fluidRow
            
            div(
              style = "margin-bottom: 22px; font-size: 89%; height: 720px; overflow: auto;",
              dataTableOutput("ctr_wells_dt")
            ) # end div
            
          ) # end tabPanel
          
        ), # end tabsetPanel
        
        fluidRow(
          
          column(
            6,
            div(selectInput(inputId = "sample_var", label = "Sampling variable:",
                            choices = c("None" = "none", "Well ID" = "file", "Drug" = "drug", "Control" = "control", "Community" = "community"),
                            selected = NULL, width = "100%"),
                style = "font-size:95%")
          ), # end column
          
          column(
            6,
            div(numericInput("max_per_group", "Events per group:", min = 1000, max = 50000, step = 1000, value = 5000),
                style = "font-size:95%")
          ) # end column
          
        ), # end fluidRow
        
        fluidRow(
          
          column(
            4,
            actionButton("plot_scatter", "Plot", width = "100%", style = "margin-bottom:20px")
          ), # end column
          
          column(
            8,
            div(textOutput("scatter_message_center"), style = "color:#c73824; margin-top:11px")
          ) # end column
          
        ) # end fluidRow
        
      ) # end column
      
    ), # end fluidRow
    
    
    fluidRow(
      
      column(
        6,
        
        tabsetPanel(
          
          tabPanel(
            title = "Scatterplot",
            br(),
            
            fluidRow(
              
              column(
                4,
                div(selectInput("scatter_var", "Point color:", choices = c("Sampling variable" = "sample_var", "Density" = "density", "None" = "none"),
                                selected = "sample_var", width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                8,
                div(selectInput("scatter_groups", "Show groups:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
            ), # end fluidRow
            
            fluidRow(
              
              column(
                4,
                div(selectInput("scatter_margins", "Margins:", choices = c("Density" = "density", "Boxplot" = "boxplot", "None" = "none"),
                                selected = "density", width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                4,
                div(selectInput("scatter_height", "Plot height:", choices = seq(300, 1500, 50), selected = 700, width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                4,
                div(selectInput("scatter_dot_size", "Point size:", choices = seq(0.1, 2, 0.1), selected = 0.5, width = "100%"),
                    style = "font-size:95%")
              ) # end column
              
            ), # end fluidRow
            
            fluidRow(
              
              column(
                6,
                div(selectInput("scatter_x", "Channel X:", choices = NULL, selected = NULL, width = "100%"),
                    style = "font-size:95%;")
              ), # end column
              
              column(
                6,
                div(selectInput("scatter_y", "Channel Y:", choices = NULL, selected = NULL, width = "100%"),
                    style = "font-size:95%")
              ) # end column
              
            ), # end fluidRow
            
            plotOutput("scatterplot", height = "700px")
            
          ) # end tabPanel
          
        ) # end tabsetPanel
        
      ), # end column
      
      column(
        6,
        
        tabsetPanel(
          
          tabPanel(
            title = "Violins",
            br(),
            
            fluidRow(
              
              column(
                4,
                div(selectInput(inputId = "violins_var", "Group by:", choices = c("Sampling variable" = "sample_var", "None" = "none"),
                                selected = "sample_var", width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                4,
                div(selectInput("violins_groups", "Show groups:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                    style = "font-size:95%")
              ), # end column
              
              column(
                4,
                div(selectInput("channels_for_violins", "Channels:", choices = NULL, selected = NULL, multiple = TRUE, width = "100%"),
                    style = "font-size:95%")
              ) # end column
              
            ), # end fluidRow
            
            plotOutput("violins_plot", height = "800px")
            
          ) # end tabPanel
          
        ) # end tabsetPanel
        
      ) # end column
      
    ) # end fluidRow
    
  ), # end tabPanel
  
  tabPanel(
    title = "Data tables", #############################################################################################
    
    tabsetPanel(
      
      tabPanel(
        title = "Drugs",
        br(),
        div(dataTableOutput("drugs_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Channels",
        br(),
        div(dataTableOutput("channels_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Channel MFIs",
        br(),
        div(dataTableOutput("mfis_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques",
        br(),
        div(dataTableOutput("cliques_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Samples edges",
        br(),
        div(dataTableOutput("samples_network_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Dispersion edges",
        br(),
        div(dataTableOutput("disp_network_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques UMAP",
        br(),
        div(dataTableOutput("umap_dt"), style = "font-size:95%")
      ), # end tabPanel
      
      tabPanel(
        title = "Cliques centroids",
        br(),
        div(dataTableOutput("centroids_dt"), style = "font-size:95%")
      ) # end tabPanel
      
    ) # end tabsetPanel
    
  ), # end tabPanel
  
  br(),
  br()
) # end ui



########################################################################################################################
########################################################################################################################
### --- Shiny SERVER --- ###############################################################################################
########################################################################################################################
########################################################################################################################

server <- function(input, output, session) { 
  
  ### set options
  options(shiny.maxRequestSize = 3000*1024^2,
          stringsAsFactors = FALSE
  ) # end options
  
  ### list of reactive values
  rvs <- reactiveValues(
    
    drugs = NULL, # data frame, drugs table
    ctr_wells = NULL, # data frame, wells from community 0
    cliques = NULL, # data frame, cliques table
    mfis = NULL, # data frame, channel MFI values in each well
    channels = NULL, # data frame, table of channel names and descriptions from FCS files
    umap = NULL, # data frame, UMAP coordinates of cliques
    centroids = NULL, # data frame, median marker expression values in cliques
    
    disp_network = NULL, # data frame, list of edges in the dispersion network
    samples_network = NULL, # data frame, list of edges in the samples network
    
    selected_wells = NULL, # vector, well ids of wells selected from UMAP and network plots
    
    sampled_event_data = NULL # data frame, down-sampled event-level data for scatterplot and violins plot
    
  ) # end rvs
  
  
  ######################################################################################################################
  ### --- Workflow tab --- #############################################################################################
  ######################################################################################################################
  
  
  ############################################################################################################
  ### --- Render UI elements for selected steps --- ##########################################################
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
  ### --- Run workflow --- ###################################################################################
  ############################################################################################################
  
  ### validate and reference file paths
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
  
  
  ### run selected workflow steps
  observeEvent(input$run_workflow, {
    
    if (!(TRUE %in% c(input$include_step1,
                      input$include_step2,
                      input$include_step3,
                      input$include_step4,
                      input$include_step5,
                      input$include_step6,
                      input$include_step7,
                      input$include_step8,
                      input$include_step9))) {
      html(id = "run_message_center",
           html = "Please check at least one workflow step to run",
           add = FALSE)
      return(NULL)
    } # end if
    
    html(id = "run_message_center",
         html = "",
         add = FALSE)
    
    
    ### RUN STEP 1
    if (input$include_step1) {
      
      validate(need(input$min_events, message = FALSE))
      
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
      
      # check if the needed files are present
      files <- list.files(paths()$out)
      
      if (!("compare_clustering.RData" %in% files)) {
        html(id = "run_message_center",
             html = "Can't find compare_clustering.RData file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      if (!("drugs_table.tsv" %in% files)) {
        html(id = "run_message_center",
             html = "Can't find drugs_table.tsv file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      if (!("Cliques.tsv" %in% files)) {
        html(id = "run_message_center",
             html = "Can't find Cliques.tsv file in output folder, please make sure that Step 8 was run",
             add = TRUE)
        return(NULL)
      } # end if
      
      html(id = "message_output", html = paste0(Sys.time(), " : Reading COMPARE files</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      channels <- strsplit(input$channels_step9, split = '[,]')[[1]]
      
      # load the compaRe RData file
      e = new.env()
      rdata_obj <- load(paste0(paths()$out, "compare_clustering.RData"), envir = e)
      rdata_obj <- e[[rdata_obj]]
      
      # drugs table
      rvs$drugs <- read.table(paste0(paths()$out, "drugs_table.tsv"), header = TRUE, sep = "\t")
      
      rvs$drugs$control <- as.character(rvs$drugs$control)
      rvs$drugs$community <- as.character(rvs$drugs$community)
      
      
      ####### samples network
      
      rvs$samples_network <- as.data.frame(cbind(as_edgelist(rdata_obj$samples_graph, names = TRUE),
                                                 as.double(E(rdata_obj$samples_graph)$weight)))
      colnames(rvs$samples_network) <- c("from","to","weight")
      
      rvs$samples_network$weight <- as.double(rvs$samples_network$weight)
      
      html(id = "message_output", html = paste0(Sys.time(), " : Assigning edge weight colors in samples network ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # option to color edges by weight for when the network is rendered
      col_fun <- colorRamp2(c(min(rvs$samples_network$weight), median(rvs$samples_network$weight), max(rvs$samples_network$weight)),
                            c("#008744", "#f1f1f1", "#a200ff"))
      rvs$samples_network$color_weight <- col_fun(rvs$samples_network$weight)
      
      # round similarity scores
      rvs$drugs$sim_vs_control <- round(rvs$drugs$sim_vs_control, 3)
      rvs$drugs$sim_vs_all <- round(rvs$drugs$sim_vs_all, 3)
      
      # calculate node degrees in samples network
      rvs$drugs$degree_samples <- NA
      for (i in 1:nrow(rvs$drugs)) {
        rvs$drugs[i,]$degree_samples <- length(which(rvs$samples_network$from == rvs$drugs[i,]$file)) + length(which(rvs$samples_network$to == rvs$drugs[i,]$file))
      } # end for
      
      
      ####### dispersion network
      
      rvs$disp_network <- as.data.frame(cbind(as_edgelist(rdata_obj$dispersion_graph, names = TRUE),
                                              as.double(E(rdata_obj$dispersion_graph)$weight)))
      colnames(rvs$disp_network) <- c("from","to","weight")
      
      rvs$disp_network$weight <- as.double(rvs$disp_network$weight)
      
      # calculate node degrees in samples network
      rvs$drugs$degree_disp <- NA
      for (i in 1:nrow(rvs$drugs)) {
        rvs$drugs[i,]$degree_disp <- length(which(rvs$disp_network$from == rvs$drugs[i,]$file)) + length(which(rvs$disp_network$to == rvs$drugs[i,]$file))
      } # end for
      
      html(id = "message_output", html = paste0(Sys.time(), " : Assigning edge weight colors in dispersion network ...</br>"), add = TRUE)
      
      # scroll window
      session$sendCustomMessage(type = "scrollCallback", 1)
      
      # option to color edges by weight for when the network is rendered
      col_fun <- colorRamp2(c(min(rvs$disp_network$weight), median(rvs$disp_network$weight), max(rvs$disp_network$weight)),
                            c("#008744", "#f1f1f1", "#a200ff"))
      rvs$disp_network$color_weight <- col_fun(rvs$disp_network$weight)
      
      
      ### read the FCS files and get the number of cells, MFI values, and channels table
      
      mfis <- matrix(data = NA, nrow = nrow(rvs$drugs), ncol = length(channels) + 1)
      colnames(mfis) <- c("well", channels)
      mfis[,"well"] <- rvs$drugs$file
      
      for (i in 1:nrow(rvs$drugs)) {
        
        well <- rvs$drugs[i,]$file
        fname <- paste0(well, ".fcs")
        
        if (fname %in% list.files(paths()$data, pattern = "fcs")) {
          
          html(id = "message_output", html = paste0(Sys.time(), " : Processing ", fname, " - file ", i, " of ", nrow(rvs$drugs), "</br>"), add = TRUE)
            
          # scroll window
          session$sendCustomMessage(type = "scrollCallback", 1)
          
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
          
        } else {
          
          html(id = "message_output", html = paste0(Sys.time(), " : Can't find FCS file for ", well), add = TRUE)
          
          # scroll window
          session$sendCustomMessage(type = "scrollCallback", 1)
          
          return(NULL)
          
        } # end else
        
      } # end for
      
      rvs$mfis <- as.data.frame(mfis)
      rvs$channels <- as.data.frame(channels_table)
      
      rownames(rvs$channels) <- rvs$channels$name
      colnames(rvs$mfis) <- c("well", rvs$channels[channels,]$desc)
      
      for (desc in rvs$channels[channels,]$desc) {
        rvs$mfis[,desc] <- as.double(rvs$mfis[,desc])
      } # end for
      
      
      ####### cliques UMAP and centroids
      
      rvs$umap <- rdata_obj$umap_
      rvs$centroids <- read.table(paste0(paths()$out, "Centroids.tsv"), sep = "\t", header = TRUE, row.names = 1)
      
      colnames(rvs$centroids) <- rvs$channels[channels,]$desc
      
      
      ####### cliques table
      
      cliques <- read.table(paste0(paths()$out, "Cliques.tsv"), sep = "\t", header = TRUE)
      
      cliques$ID <- paste0("C", cliques$ID) # to match cliques ID in umap and centroids tables
      cliques$length <- NA # number of wells in clique
      cliques$total_live_cells <- NA # total events in clique
      
      for (i in 1:nrow(cliques)) {
        cliques[i,]$length <- length(strsplit(cliques[i,]$File, split = ",")[[1]])
        cliques[i,]$total_live_cells <- sum(as.numeric(strsplit(cliques[i,]$live_cels, split = ",")[[1]]))
        
      } # end for 
      
      rvs$drugs$clique_ids <- NA
      
      for (i in 1:nrow(rvs$drugs)) {
        
        for (c in 1:nrow(cliques)) {
          
          if (rvs$drugs$file[i] %in% strsplit(cliques[c,]$File, split = ",")[[1]]) {
            
            if (is.na(rvs$drugs[i,]$clique_ids))
              rvs$drugs[i,]$clique_ids <- cliques[c,]$ID
            else
              rvs$drugs[i,]$clique_ids <- paste(c(rvs$drugs[i,]$clique_ids, cliques[c,]$ID), collapse = ",")
            
          } # end if
          
        } # end for
        
        
      } # end for
      
      rvs$cliques <- cliques
      rvs$ctr_wells <- rvs$drugs[which(rvs$drugs$community == "0"),]
      
      ####### update ui controls
      
      updateSelectInput(session, "umap_x", choices = colnames(rvs$umap))
      updateSelectInput(session, "umap_y", choices = colnames(rvs$umap), selected = colnames(rvs$umap)[2])
      updateSelectInput(session, "umap_channel", choices = as.character(rvs$channels[channels,]$desc))
      
      updateSelectInput(session, "samples_node_color_channel", choices = as.character(rvs$channels[channels,]$desc))
      updateSelectInput(session, "disp_node_color_channel", choices = as.character(rvs$channels[channels,]$desc))
      
      min_degree <- min(rvs$drugs$degree_samples)
      max_degree <- max(rvs$drugs$degree_samples)
      
      updateSliderInput(session, "node_degree_range", min = min(rvs$drugs$degree_samples), max = max(rvs$drugs$degree_samples),
                        value = c(min(rvs$drugs$degree_samples), max(rvs$drugs$degree_samples)))
      
      min_weight <- round(min(rvs$samples_network$weight), 1) - 1
      max_weight <- round(max(rvs$samples_network$weight), 1) + 1
      
      updateSliderInput(session, "edge_weight_range", min = min_weight, max = max_weight, value = c(min_weight, max_weight))
      
      min_score <- round(min(rvs$drugs$sim_vs_control), 1) - 1
      max_score <- round(max(rvs$drugs$sim_vs_control), 1) + 2
      
      updateSliderInput(session, "sim_vs_ctr_range", min = min_score, max = max_score, value = c(min_score, max_score))
      
      min_score <- round(min(rvs$drugs$sim_vs_all), 1) - 1
      max_score <- round(max(rvs$drugs$sim_vs_all), 1) + 1
      
      updateSliderInput(session, "sim_vs_all_range", min = min_score, max = max_score, value = c(min_score, max_score))
      
      names(channels) <- rvs$channels[channels,]$desc
      
      updateSelectInput(session, "scatter_x", choices = channels)
      updateSelectInput(session, "scatter_y", choices = channels, selected = channels[2])
      updateSelectInput(session, "channels_for_violins", choices = channels, selected = channels)
      
      
      html(id = "message_output", html = paste0(Sys.time(), " : Step 9 done! Proceed to Plots tab.</br>"), add = TRUE)
      session$sendCustomMessage(type = "scrollCallback", 1)
      
    } # end if
    
  }) # end observeEvent
  
  
  
  ######################################################################################################################
  ### --- Data tables tab --- ##########################################################################################
  ######################################################################################################################
  
  ### render the drugs table
  output$drugs_dt <- renderDT(
    rvs$drugs[,c("file","drug","concentration","live_cells","control","community","sim_vs_control","sim_vs_all","clique_ids")],
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$drugs_dt
  
  ### render the channels table
  output$channels_dt <- renderDT(
    rvs$channels,
    selection = "none",
    rownames = TRUE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$channels_dt
  
  ### render the channel MFIs table
  output$mfis_dt <- renderDT(
    rvs$mfis,
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$mfis_dt
  
  ### render the cliques table
  output$cliques_dt <- renderDT(
    rvs$cliques,
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$cliques_dt
  
  ### render the samples network table
  output$samples_network_dt <- renderDT(
    rvs$samples_network,
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$samples_network_dt
  
  ### render the samples network table
  output$disp_network_dt <- renderDT(
    rvs$disp_network,
    filter = "top",
    selection = "none",
    rownames = FALSE,
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$disp_network_dt
  
  ### render the cliques UMAP table
  output$umap_dt <- renderDT(
    rvs$umap,
    selection = "none",
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$umap_dt
  
  ### render the cliques centroids table
  output$centroids_dt <- renderDT(
    rvs$centroids,
    selection = "none",
    options = list(pageLength = 25, lengthMenu = c(10, 25, 50, 100, 500, 1000), searchHighlight = TRUE)
  ) # end output$umap_dt
  
  
  
  ######################################################################################################################
  ### --- Plots tab --- ###############################################################################################
  ######################################################################################################################
  
  
  ### render the drugs table
  output$ctr_wells_dt <- renderDT(
    rvs$ctr_wells[,c("file","drug","concentration","control","live_cells")],
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
  
  
  ############################################################################################################
  ### --- make and render the UMAP scatterplot --- ###########################################################
  ############################################################################################################
  
  html(id = "umap_plot_hover_name",
       html = paste("<strong>Clique:</strong>"),
       add = FALSE)
  html(id = "umap_plot_hover_centr",
       html = paste("<strong>MFI:</strong>"),
       add = FALSE)
  html(id = "umap_plot_hover_wells",
       html = paste("<strong>Wells:</strong>"),
       add = FALSE)
  html(id = "umap_plot_hover_drugs",
       html = paste("<strong>Drugs:</strong>"),
       add = FALSE)
  
  ### data for the umap scatterplot
  umap_plot_data <- reactive({
    
    validate(need(rvs$umap, message = FALSE))
    validate(need(rvs$centroids, message = FALSE))
    validate(need(rvs$cliques, message = FALSE))
    validate(need(input$umap_x, message = FALSE))
    validate(need(input$umap_y, message = FALSE))
    validate(need(input$umap_channel, message = FALSE))
    
    d <- data.frame(
      x = rvs$umap[,input$umap_x],
      y = rvs$umap[,input$umap_y],
      centr = rvs$centroids[,input$umap_channel],
      name = rownames(rvs$umap),
      wells = c("", rvs$cliques$File),
      drugs = c("", rvs$cliques$Clique),
      comm = c("Control", rvs$cliques$Community)
    ) # end data.frame
    
    if (input$umap_dot_size == "total_live_cells") {
      d$size <- c(max(rvs$cliques$total_live_cells), rvs$cliques$total_live_cells)
      size_title <- "Total number of events"
    } else if (input$umap_dot_size == "length") {
      d$size <- c(max(rvs$cliques$length), rvs$cliques$length)
      size_title <- "Number of wells"
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
    
    border_strokes <- rep(0, nrow(d))
    border_strokes[which(d$name == "Control")] <- 2
    
    border_cols <- rep("#ffffff", nrow(d))
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
        guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 16, angle = 45, hjust = 1)),
               size = guide_legend(order = 2, title = umap_plot_data()$size_title, label.theme = element_text(size = 16))) +
        theme_light() +
        theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.title = element_text(size = 20), axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
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
        guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 16, angle = 45, hjust = 1)),
               size = guide_legend(order = 2, title = umap_plot_data()$size_title, label.theme = element_text(size = 16))) +
        theme_light() +
        theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.title = element_text(size = 20), axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        scale_fill_gradientn(colours = c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
                             limits = c(min(d$centr), max(d$centr)),
                             breaks = round(seq(min(d$centr) + 0.01, max(d$centr) - 0.01, length.out = 3), 2))
      
    } else {
      
      p <- ggplot(d, aes(x = x, y = y)) +
        geom_point(aes(fill = centr), size = 10, shape = node_shapes, alpha = 0.8, color = border_cols, stroke = border_strokes) +
        labs(x = umap_plot_data()$x_lab, y = umap_plot_data()$y_lab) +
        guides(fill = guide_colorbar(order = 1, title = umap_plot_data()$channel, label.theme = element_text(size = 16, angle = 45, hjust = 1))) +
        theme_light() +
        theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.title = element_text(size = 20), axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        scale_fill_gradientn(colours = c("#0571B0","#92C5DE","yellow","#F4A582","#CA0020"),
                             limits = c(min(d$centr), max(d$centr)),
                             breaks = round(seq(min(d$centr) + 0.01, max(d$centr) - 0.01, length.out = 3), 2))
       
    } # end else
    
    if (input$umap_labels == "Clique ID") {
      
      p <- p + geom_text_repel(box.padding = unit(0.3, "lines"), label = d$name,
                               size = 5.5, show.legend = FALSE, fontface = 1)
      
    } else if (input$umap_labels == "Community") {
      
      color_palette <- GetColors(length(unique(d$comm)), scheme = "jet")
      
      p <- p + geom_text_repel(box.padding = unit(0.6, "lines"), label = d$comm,
                               size = 7, show.legend = FALSE, fontface = 2, color = color_palette[as.factor(d$comm)])
      
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
    
    if (nrow(point) > 0) {
      wells <- strsplit(point$wells, split = ",")[[1]]
      rvs$selected_wells <- union(rvs$selected_wells, wells)
    } # end if
    
  }) # end observeEvent
  
  
  ### subset of annotation with wells clicked in line plot and plate heatmaps
  scatter_wells <- reactive({
    if (is.null(rvs$selected_wells))
      return(NULL)
    else {
      rows_to_select <- which(rvs$drugs$file %in% rvs$selected_wells & !(rvs$drugs$file %in% rvs$ctr_wells$file))
      validate(need(rows_to_select, message = FALSE))
      return(rvs$drugs[rows_to_select,])
    } # end else
  }) # end scatter_wells
  
  
  ### observer to render the datatable with wells selected from line plot and plate heatmaps
  observe({
    
    if (is.null(scatter_wells())) {
      output$scatter_wells_dt <- renderDT(NULL)
    } else {
      output$scatter_wells_dt <- renderDT(
        scatter_wells()[,c("file","drug","concentration","community","clique_ids","live_cells")],
        rownames = FALSE,
        options = list(pageLength = 500, lengthMenu = c(10,25,50,100,500,1000), searchHighlight = TRUE)
      ) # end output$scatter_wells_dt
    } # end else
    
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
  observeEvent(input$scatter_wells_clear_table, {
    
    isolate(visNetworkProxy("disp_network") %>% visGetSelectedNodes())
    isolate(visNetworkProxy("samples_network") %>% visGetSelectedNodes())
    
    rvs$selected_wells <- NULL
  }) # end observeEvent
  
  
  ############################################################################################################
  ### --- make and render the dispersion network --- #########################################################
  ############################################################################################################
  
  ### samples graph igraph object
  disp_igraph <- reactive({
    
    validate(need(rvs$drugs, message = FALSE))
    
    nodes <- data.frame(
      id = rvs$drugs$file,
      color.border = "#999999",
      borderWidth = 1,
      color.highlight.border = "#d62d20",
      degree = rvs$drugs$degree_disp,
      value = rvs$drugs$live_cells,
      drug = rvs$drugs$drug,
      drug_label = paste(rvs$drugs$drug, rvs$drugs$concentration, sep = "_"),
      community = rvs$drugs$community,
      sim_vs_ctr = rvs$drugs$sim_vs_control,
      sim_vs_all = rvs$drugs$sim_vs_all
    ) # end data.frame 
    
    # assign node color based on selected metric
    if (input$disp_node_color_metric == "none") {
      
      nodes$color.background <- "#f1f1f1"
      nodes$color.highlight.background <- "#f1f1f1"
      
    } else {
      
      if (input$disp_node_color_metric == "degree")
        col_values <- rvs$drugs$degree_disp
      else if (input$disp_node_color_metric == "cells")
        col_values <- rvs$drugs$live_cells
      else if (input$disp_node_color_metric == "mfi")
        col_values <- rvs$mfis[,input$disp_node_color_channel]
      else if (input$disp_node_color_metric == "sim_vs_ctr")
        col_values <- rvs$drugs$sim_vs_control
      else if (input$disp_node_color_metric == "sim_vs_all")
        col_values <- rvs$drugs$sim_vs_all
      
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
                                     degree = length(which(rvs$disp_network$from == "Control")) + length(which(rvs$disp_network$to == "Control")),
                                     value = mean(rvs$drugs[which(rvs$drugs$degree_disp == 0),]$live_cells),
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
        col_values <- c(rvs$drugs$degree_disp, ctr_value)
        
      } else if (input$disp_node_color_metric == "cells") {
        
        ctr_value <- nodes[which(nodes$id == "Control"),]$value
        col_values <- c(rvs$drugs$live_cells, ctr_value)
        
      } else if (input$disp_node_color_metric == "mfi") {
        
        ctr_value <- median(rvs$mfis[which(rvs$drugs$degree_disp == 0),input$disp_node_color_channel])
        col_values <- c(rvs$mfis[,input$disp_node_color_channel], ctr_value)
        
      } else if (input$disp_node_color_metric == "sim_vs_ctr") {
        
        ctr_value <- mean(rvs$drugs[which(rvs$drugs$degree_disp == 0),]$sim_vs_control)
        col_values <- c(rvs$drugs$sim_vs_control, ctr_value)
        
      } else if (input$disp_node_color_metric == "sim_vs_all") {
        
        ctr_value <- mean(rvs$drugs[which(rvs$drugs$degree_disp == 0),]$sim_vs_all)
        col_values <- c(rvs$drugs$sim_vs_all, ctr_value)
        
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
      from = rvs$disp_network$from,
      to = rvs$disp_network$to,
      weight = rvs$disp_network$weight
    ) # end data.frame
    
    # assign edge colors
    if (input$disp_edge_color_metric == "weight") {
      edges$color <- rvs$disp_network$color_weight
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
    updateSelectInput(session, "disp_select_by_comm", choices = V(g)$community[order(V(g)$community)])
    
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
        visOptions(highlightNearest = input$disp_highlight_nearest, autoResize = TRUE, ) %>%
        visInteraction(multiselect = TRUE)
      
    } # end else if
    
  }) # end output$disp_network
  
  
  ### observer to add selected nodes form the dispersion network into the scatter well table
  observe({
    
    input$disp_network_add_selected
    
    visNetworkProxy("disp_network") %>% visGetSelectedNodes()
    
    validate(need(input$disp_network_selectedNodes, message = FALSE))
    
    rvs$selected_wells <- union(isolate(rvs$selected_wells), setdiff(input$disp_network_selectedNodes, "Control"))
    
  }) # end observe
  
  
  ############################################################################################################
  ### --- observers to select nodes in the dispersion network --- ############################################
  ############################################################################################################
  
  ### observer to select nodes by id
  observeEvent(input$disp_select_by_id, ignoreNULL = FALSE, {
    
    if (is.null(input$disp_select_by_id))
      visNetworkProxy("disp_network") %>% visUnselectAll()
    else
      visNetworkProxy("disp_network") %>% visSelectNodes(id = input$disp_select_by_id)
    
  }) # end observeEvent
  
  
  ### observer to select nodes by drug
  observeEvent(input$disp_select_by_drug, ignoreNULL = FALSE, {
    
    if (is.null(input$disp_select_by_drug))
      visNetworkProxy("disp_network") %>% visUnselectAll()
    else {
      g <- disp_igraph()
      visNetworkProxy("disp_network") %>%
        visSelectNodes(id = V(g)$name[which(V(g)$drug %in% input$disp_select_by_drug)])
    } # end else
    
  }) # end observeEvent
  
  
  ### observer to select nodes by community
  observeEvent(input$disp_select_by_comm, ignoreNULL = FALSE, {
    
    if (is.null(input$disp_select_by_comm))
      visNetworkProxy("disp_network") %>% visUnselectAll()
    else {
      g <- disp_igraph()
      visNetworkProxy("disp_network") %>%
        visSelectNodes(id = V(g)$name[which(V(g)$community %in% input$disp_select_by_comm)]) 
    } # end else
    
  }) # end observeEvent
  
  
  ############################################################################################################
  ### --- make and render the samples network --- ############################################################
  ############################################################################################################
  
  ### samples graph igraph object
  samples_igraph <- reactive({
    
    validate(need(rvs$drugs, message = FALSE))
    
    nodes <- data.frame(
      id = rvs$drugs$file,
      color.border = "#999999",
      borderWidth = 1,
      color.highlight.border = "#d62d20",
      value = rvs$drugs$live_cells,
      degree = rvs$drugs$degree_samples,
      sim_vs_ctr = rvs$drugs$sim_vs_control,
      sim_vs_all = rvs$drugs$sim_vs_all,
      community = rvs$drugs$community,
      control = rvs$drugs$control,
      drug = rvs$drugs$drug,
      clique_ids = rvs$drugs$clique_ids
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
        col_values <- rvs$drugs$degree_samples
      else if (input$samples_node_color_metric == "cells")
        col_values <- rvs$drugs$live_cells
      else if (input$samples_node_color_metric == "mfi")
        col_values <- rvs$mfis[,input$samples_node_color_channel]
      else if (input$samples_node_color_metric == "sim_vs_ctr")
        col_values <- rvs$drugs$sim_vs_control
      else if (input$samples_node_color_metric == "sim_vs_all")
        col_values <- rvs$drugs$sim_vs_all
      
      col_fun <- colorRamp2(c(min(col_values), median(col_values), max(col_values)),
                            c("#0057e7", "#fdf498", "#d62d20"))
      nodes$color.background <- col_fun(col_values)
      nodes$color.highlight.background <- col_fun(col_values)
      
    } # end else
    
    # assign node label based on user input
    if (input$samples_node_label == "well") {
      nodes$label <- rvs$drugs$file
    } else if (input$samples_node_label == "drug") {
      
      if ("concentration" %in% colnames(rvs$drugs))
        nodes$label <- paste(rvs$drugs$drug, rvs$drugs$concentration, sep = "_")
      else
        nodes$label <- rvs$drugs$drug
      
    } else if (input$samples_node_label == "comm") {
      nodes$label <- rvs$drugs$community
    } # end esle if
    
    # make edges data frame
    edges <- data.frame(
      from = rvs$samples_network$from,
      to = rvs$samples_network$to,
      weight = rvs$samples_network$weight
    ) # end data.frame
    
    # assign edge colors
    if (input$samples_edge_color_metric == "weight") {
      edges$color <- rvs$samples_network$color_weight
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
    updateSelectInput(session, "samples_select_by_comm", choices = V(g)$community[order(V(g)$community)])
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
  
  
  ### observer to add selected nodes form the samples network into the scatter well table
  observe({
    
    input$samples_network_add_selected
    
    visNetworkProxy("samples_network") %>% visGetSelectedNodes()
    
    validate(need(input$samples_network_selectedNodes, message = FALSE))
    
    rvs$selected_wells <- union(isolate(rvs$selected_wells), input$samples_network_selectedNodes)
    
  }) # end observe
  
  
  ############################################################################################################
  ### --- observers to select nodes in the samples network --- ###############################################
  ############################################################################################################
  
  ### observer to select nodes by id
  observeEvent(input$samples_select_by_id, ignoreNULL = FALSE, {
    
    if (is.null(input$samples_select_by_id))
      visNetworkProxy("samples_network") %>% visUnselectAll()
    else
      visNetworkProxy("samples_network") %>% visSelectNodes(id = input$samples_select_by_id)
    
  }) # end observeEvent
  
  
  ### observer to select nodes by drug
  observeEvent(input$samples_select_by_drug, ignoreNULL = FALSE, {
    
    if (is.null(input$samples_select_by_drug))
      visNetworkProxy("samples_network") %>% visUnselectAll()
    else {
      g <- samples_igraph()
      visNetworkProxy("samples_network") %>%
        visSelectNodes(id = V(g)$name[which(V(g)$drug %in% input$samples_select_by_drug)])
    } # end else
    
  }) # end observeEvent
  
  
  ### observer to select nodes by community
  observeEvent(input$samples_select_by_comm, ignoreNULL = FALSE, {
    
    if (is.null(input$samples_select_by_comm))
      visNetworkProxy("samples_network") %>% visUnselectAll()
    else {
      g <- samples_igraph()
      visNetworkProxy("samples_network") %>%
        visSelectNodes(id = V(g)$name[which(V(g)$community %in% input$samples_select_by_comm)]) 
    } # end else
    
  }) # end observeEvent
  
  
  ### observer to select nodes by clique
  observeEvent(input$samples_select_by_clique, ignoreNULL = FALSE, {
    
    if (is.null(input$samples_select_by_clique))
      visNetworkProxy("samples_network") %>% visUnselectAll()
    else {
      g <- samples_igraph()
      
      ids <- NULL
      for (c in input$samples_select_by_clique) {
        ids <- union(ids, strsplit(rvs$cliques[which(rvs$cliques$ID == c),]$File, split = ",")[[1]])
      } # end for
      
      visNetworkProxy("samples_network") %>%
        visSelectNodes(id = ids) 
    } # end else
    
  }) # end observeEvent
  
  
  ############################################################################################################
  ### --- scatterplot --- ####################################################################################
  ############################################################################################################
  
  
  ### collect data for scatterplot and violins plot
  observeEvent(input$plot_scatter, {
    
    rvs$sampled_event_data <- sampled_event_data()
    
  }) # end observeEvent
  
  
  ### event expression data
  event_data <- reactive({
    
    validate(need(scatter_wells(), message = FALSE))
    
    wells <- rbind(scatter_wells()[input$scatter_wells_dt_rows_selected,], rvs$ctr_wells[input$ctr_wells_dt_rows_selected,])
    
    channel_names <- strsplit(input$channels_step9, split = '[,]')[[1]]
    
    html(id = "scatter_message_center", html = "Reading FCS files ...", add = FALSE)
    
    files <- paste0(paths()$data, wells$file, ".fcs")
    
    # read fcs files into FlowSet
    fset_obj <- read.flowSet(files, transformation = FALSE)
    
    html(id = "scatter_message_center", html = "Applying transformation ...", add = FALSE)
    
    # apply transformation
    fset_obj <- transform(fset_obj, transformList(channel_names, logicleTransform(w = 0.5, t = 262144, m = 4.5)))
    
    vars <- c("file","drug","control","community")
    
    # extract into data frame and record the start and end row positions of each well in the exprs_data
    html(id = "scatter_message_center", html = "Extracting values into data frame ...", add = FALSE)
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
    
    # center and scale
    html(id = "scatter_message_center", html = "Centering and scaling ...", add = FALSE)
    exprs_data  <- scale(exprs_data, center = TRUE, scale = TRUE)
    
    # attach color columns for variables
    html(id = "scatter_message_center", html = "Generating color palettes for variables ...", add = FALSE)
    var_color_cols <- NULL
    for (c in 1:ncol(var_cols)) {
      color_palette <- GetColors(length(unique(var_cols[,c])), scheme = "smooth rainbow", start = 0.2, end = 0.9, bias = 1.6)
      var_color_cols <- cbind(var_color_cols, color_palette[as.factor(var_cols[,c])])
    } # end for
    
    exprs_data <- as.data.frame(cbind(exprs_data, var_cols, var_color_cols))
    
    colnames(exprs_data) <- c(channel_names, vars, paste0(vars, "_color"))
    
    for (i in 1:length(channel_names)) {
      exprs_data[,i] <- as.double(exprs_data[,i])
    } # end for
    
    #write.table(exprs_data, "event_exprs_data.txt", sep = "\t", quote = FALSE)
    
    html(id = "scatter_message_center", html = "", add = FALSE)
    
    return(exprs_data)
    
  }) # end reactive
  
  
  ### event expression data after downsampling
  sampled_event_data <- reactive({
    
    validate(need(event_data(), message = FALSE))
    
    html(id = "scatter_message_center", html = "Sampling ...", add = FALSE)
    
    d <- event_data()
    
    sv <- isolate(input$sample_var)
    m <- isolate(input$max_per_group)
    
    if (sv == "none") {
      
      if (m >= nrow(d)) {
        d_sampled <- d
      } else {
        d_sampled <- d[sample(nrow(d), m),]
        d_sampled <- d_sampled[order(as.numeric(rownames(d_sampled))),]
      } # end else
      
    } else {
      
      var_groups <- unique(d[,sv])
      
      d_sampled <- NULL
      
      for (group in var_groups) {
        
        d_group <- d[which(d[,sv] == group),]
        
        if (m >= nrow(d_group)) {
          d_sampled <- rbind(d_sampled, d_group)
        } else {
          d_group_sampled <- d_group[sample(nrow(d_group), m),]
          d_group_sampled <- d_group_sampled[order(as.numeric(rownames(d_group_sampled))),]
          d_sampled <- rbind(d_sampled, d_group_sampled)
        } # end else
        
      } # end for 
      
    } # end else
    
    html(id = "scatter_message_center", html = "", add = FALSE)
    
    #write.table(d_sampled, "event_exprs_data_sampled.txt", sep = "\t", quote = FALSE)
    
    return(d_sampled)
    
  }) # end sampled_event_data()
  
  
  ### change scatterplot group depending on the scatterplot variable
  observe({
    
    input$scatter_var
    
    if (input$scatter_var == "sample_var") {
      updateSelectInput(session, "scatter_groups",
                        choices = unique(rvs$sampled_event_data[,isolate(input$sample_var)]),
                        selected = unique(rvs$sampled_event_data[,isolate(input$sample_var)]))
    } else if (input$scatter_var == "density") {
      updateSelectInput(session, "scatter_groups", choices = c("No groups to show" = "density"), selected = "density")
    } else if (input$scatter_var == "none") {
      updateSelectInput(session, "scatter_groups", choices = c("No groups to show" = "none"), selected = "none")
    } # end else
    
  }) # end observeEvent
  
  
  ### data for scatterplot
  scatterplot_data <- reactive({
    
    validate(need(input$scatter_groups, message = FALSE))
    
    if (input$scatter_x == input$scatter_y)
      return(NULL)
    
    d <- rvs$sampled_event_data
    
    if ("none" %in% input$scatter_groups) {
      
      d <- d[,c(input$scatter_x, input$scatter_y)]
      
      return(d)
      
    } else if ("density" %in% input$scatter_groups) {
      
      d <- d[,c(input$scatter_x, input$scatter_y)]
      
      # get density in 2d for the x and y points
      d$density <- get_density(d[,1], d[,2], h = c(1,1), n = 100)
      
      return(d)
      
    } else {
      
      scatter_var <- isolate(input$sample_var)
      
      if (!all(input$scatter_groups %in% unique(d[,scatter_var])))
        return(NULL)
      
      d <- d[which(d[,scatter_var] %in% input$scatter_groups),]
      
      d_scatter <- d[,c(input$scatter_x, input$scatter_y)]
      
      # get density in 2d for the x and y points
      d_scatter$density <- get_density(d_scatter[,1], d_scatter[,2], h = c(1,1), n = 100)
      
      # attach variable and variable color to the scatterplot data
      d_scatter$var_name <- d[,scatter_var]
      d_scatter$var_color <- d[,paste0(scatter_var, "_color")]
      
      return(d_scatter)
      
    } # end else
    
  }) # end scatterplot_data
  
  
  ### return the selected height of the static scatterplot
  scatterplot_height <- function() {
    return(as.integer(input$scatter_height))
  } # end scatterplot_height
  
  ### function to make the static scatterplot
  output$scatterplot <- renderPlot(height = scatterplot_height, {
    
    validate(need(scatterplot_data(), message = FALSE))
    
    d <- scatterplot_data()
    
    x_lab <- colnames(d)[1]
    y_lab <- colnames(d)[2]
    
    colnames(d)[1] <- "x"
    colnames(d)[2] <- "y"
    
    scatter_var <- isolate(input$scatter_var)
    
    if (scatter_var == "none") { # no point colors
      
      p <- ggplot(d, aes(x = x, y = y)) +
        geom_point(size = as.double(input$scatter_dot_size), alpha = 0.7) +
        labs(x = x_lab, y = y_lab) +
        theme_light() +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20))
      
      if (input$scatter_margins == "none")
        p
      else
        ggMarginal(p, type = input$scatter_margins, size = 5)
      
    } else if (scatter_var == "density") { # color by density
      
      p <- ggplot(d, aes(x = x, y = y, color = density)) +
        geom_point(size = as.double(input$scatter_dot_size), alpha = 0.8) +
        labs(x = x_lab, y = y_lab) +
        theme_light() +
        theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.text = element_text(size = 14, angle = 45, hjust = 1), legend.title = element_text(size = 14),
              axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        scale_color_viridis_c()
      
      if (input$scatter_margins == "none")
        p
      else
        ggMarginal(p, type = input$scatter_margins, size = 5)
      
    } else { # color by all other variables
      
      p <- ggplot(d, aes(x = x, y = y, color = var_name)) +
        geom_point(size = as.double(input$scatter_dot_size), alpha = 0.8) +
        labs(x = x_lab, y = y_lab) +
        theme_light() +
        theme(legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.title = element_blank(), legend.text = element_text(size = 14),
              axis.text = element_text(size = 18), axis.title = element_text(size = 20)) +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        scale_color_manual(values = setNames(unique(d$var_color), unique(d$var_name)))
      
      if (input$scatter_margins == "none")
        p
      else
        ggMarginal(p, type = input$scatter_margins, size = 5, groupColour = TRUE, groupFill = TRUE)
      
    } # end else
    
  }) # end output$scatter_static
  
  
  ### change violins group depending on the violins variable
  observe({
    
    input$violins_var
    
    if (input$violins_var == "sample_var") {
      updateSelectInput(session, "violins_groups",
                        choices = unique(rvs$sampled_event_data[,isolate(input$sample_var)]),
                        selected = unique(rvs$sampled_event_data[,isolate(input$sample_var)]))
    } else if (input$violins_var == "none") {
      updateSelectInput(session, "violins_groups", choices = c("No groups to show" = "none"), selected = "none")
    } # end else
    
  }) # end observeEvent
  
  
  ### data for the channel intensities plot
  violins_data <- reactive({
    
    validate(need(input$channels_for_violins, message = FALSE))
    validate(need(input$violins_groups, message = FALSE))
    
    d <- rvs$sampled_event_data
    
    if ("none" %in% input$violins_groups) {
      
      data_long <- NULL
      for (c in input$channels_for_violins)
        data_long <- rbind(data_long, cbind(d[,c], rep(c, nrow(d))))
      
      data_long <- as.data.frame(data_long)
      
      colnames(data_long) <- c("intensity", "channel")
      
      data_long$intensity <- as.double(data_long$intensity)
      data_long$channel <- as.factor(data_long$channel)
      
    } else {
      
      d <- d[which(d[,isolate(input$sample_var)] %in% input$violins_groups),]
      
      data_long <- NULL
      for (c in input$channels_for_violins)
        data_long <- rbind(data_long, cbind(d[,c], rep(c, nrow(d)),
                                            d[,isolate(input$sample_var)],
                                            d[,paste0(isolate(input$sample_var), "_color")]))
      
      data_long <- as.data.frame(data_long)
      
      colnames(data_long) <- c("intensity", "channel", "var_name", "var_color")
      
      data_long$intensity <- as.double(data_long$intensity)
      data_long$channel <- as.factor(data_long$channel)
      data_long$var_name <- as.factor(data_long$var_name)
      data_long$var_color <- as.character(data_long$var_color)
      
    } # end else
    
    return(data_long)
    
  }) # end violins_data()
  
  
  ### render the channel intensities plot
  output$violins_plot <- renderPlot({
    
    validate(need(violins_data(), message = FALSE))
    
    d <- violins_data()
    
    violins_var <- isolate(input$violins_var)
    
    if (violins_var == "none") {
      
      ggplot(d, aes(x = channel, y = intensity)) +
        geom_violin(trim = FALSE, fill = "#f1f1f1") +
        geom_boxplot(width = 0.1, fill = "#f1f1f1") +
        labs(x = NULL) +
        theme_light() +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        scale_x_discrete(limits = as.factor(input$channels_for_violins))
      
    } else {
      
      ggplot(d, aes(x = channel, y = intensity, fill = var_name)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1, position = position_dodge(width = 0.9)) +
        labs(x = NULL) +
        theme_light() +
        theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.direction = "horizontal", legend.position = c("bottom"), legend.justification = c(-0.03,0),
              legend.text = element_text(size = 14), legend.title = element_blank()) +
        scale_fill_manual(values = setNames(unique(d$var_color), unique(d$var_name))) +
        scale_x_discrete(limits = as.factor(input$channels_for_violins))
      
    } # end else
    
  }) # end output$violins_plot
  
  
} # end server

### make Shiny app
shinyApp(ui = ui, server = server)
