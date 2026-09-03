library(tidyverse)
library(readxl)
library(dr4pl)
library(shiny)

# ===========================================
# UI 
# ===========================================

ui <- fluidPage(
  titlePanel("Elisa data plot and extrapolation"),
  sidebarLayout(
        # File load
        sidebarPanel(
          fileInput("file1",
                    "Choose xlsx file",
                    accept = ".xlsx"),
          
        # Cell ranges
          textInput("measurement", 
                    label = "Cell range of measurements",
                    placeholder = "e.g. A1:F97"),
        
          textInput("metadata", 
                    label = "Cell range of sample data",
                    placeholder = "e.g. H1:I20"),
        
        #Sheet name:
          textInput(inputId = "sheetname", 
                    label = "Sheetname with ELISA summary data",
                    placeholder = "e.g. Result summary"),
        
        # Slider to set extrapolation factor
          sliderInput("extrapolation", 
                      label = "Extrapolation factor",
                      min = 0,
                      max = 0.5, 
                      value = 0.1,
                      step = 0.01),
        
        # Button to export csv with calculated values.
          downloadButton("export_csv", 
                         label = "Export calculated data as CSV")
          ),
        
  
        # Implementation of tabs for all data
        mainPanel(tabsetPanel(
                    tabPanel(
                      "Standard data",
                      h4("Raw input data"),
                      tableOutput("table_raw"),
                      
                      h4("table_metadata"),
                      tableOutput("table_metadata"),
                      
                      h4("Matched standard data"),
                      tableOutput("table_standard")
                    ),
                    
                    tabPanel(
                      "Standard curve",
                      plotOutput("standard_curve"),
                      br(),
                      h4("Model parameters"),
                      verbatimTextOutput("model")
                      ),
                    
                    tabPanel(
                      "Calculated values",
                      tableOutput("table_calculated")
                  )
                )
              )
        )
  )


# ========================================================================
# Server
# ========================================================================

server <- function(input, output, session) {
  ############################################
  # Helper functions
  ############################################
  
  invert_4pl <- function(y, lower, upper, ic50, slope){
    
    value <- ((lower - upper) / (y - upper)) - 1
    
    if(!is.finite(value) || value <= 0) {
      return(NA_real_)
    }
    
    ic50 * value^(1/slope)
  }
  
  ###########################################
  # READ RAW DATA
  ###########################################
  
  raw_data <- reactive({
    req(
      input$file1,
      input$measurement,
      input$sheetname
    )
    
    validate(
      need(input$measurement != "",
           "Please enter the measurement range"),
      need(input$sheetname != "", 
           "Please enter the sheet name")
    )
    
    read_xlsx(input$file1$datapath,
              range = input$measurement,
              sheet = input$sheetname,
              col_names = c("well", "group", "sample", 
                            "absorbance450", "absorbande650",
                            "basic_calculation")
    )
  })
  
  ###########################################
  # READ CONCENTRATION METADATA
  ###########################################
  
  concentration_data <- reactive({
    req(input$file1,
        input$metadata,
        input$sheetname)
    
    validate(
      need(
        input$metadata != "",
        "Please enter the metadata range"
      )
    )
    
    read_xlsx(
      input$file1$datapath,
      range = input$metadata,
      sheet = input$sheetname,
      col_names = c("sample", 
                    "concentration")
    ) %>% 
      mutate(concentration = suppressWarnings(as.numeric(concentration)))
  })
  
  ###########################################
  # STANDARD DATA
  ###########################################
  
  standard_data <- reactive({
    standards <- raw_data() %>% 
      filter(str_detect(sample, "Std")) %>% 
      left_join(
        concentration_data(),
        by = "sample"
      ) %>% 
      drop_na(concentration, basic_calculation)
    
    validate(
      need(nrow(standards) > 0,
           "No valid standard data found")
    )
    
    standards
  })
  
  ###########################################
  # STANDARD MEANS
  ###########################################
  
  standard_means <- reactive({
    
  standard_data() %>% 
      group_by(concentration) %>% 
      summarise(mean = mean(basic_calculation),
                .groups = "drop") %>% 
      arrange(concentration)
  })
  
  ###########################################
  # FIT 4PL MODEL
  ###########################################
  
  model <- reactive({
    df <- standard_data()
    
    validate(
      need(
        n_distinct(df$concentration) >= 4,
        "At least four different standard concentrations are recommended"
      )
    )
    
    dr4pl(
      basic_calculation ~ concentration,
      data = df,
      method.init = "logistic"
    )
  })
  
  ###########################################
  # MODEL PARAMETERS
  ###########################################
  
  model_parameters <- reactive({
    params <- model()$parameters
    
    list(upper = params[1],
         ic50 = params [2],
         slope = params[3],
         lower = params[4])
  })
  
  ###########################################
  # EXTRAPOLATION MODELS
  ###########################################
 
   extrapolation_models <- reactive({
    means <- standard_means()
    
    validate(
      need(nrow(means) >= 2,
           "At least two standard concentrations are required")
    )
    
    lowest_points <- means %>% 
      slice_head(n = 2)
    
    highest_points <- means %>% 
      slice_tail(n = 2)
    
    list(
      lowest_concentration_model = 
        lm(mean ~ concentration, data = lowest_points),
      
      highest_concentration_model =
        lm(mean ~ concentration, data = highest_points),
      
      lowest_points = lowest_points,
      highest_points = highest_points
    )
  })
  
  ###########################################
  # EXTRAPOLATION LIMITS
  ###########################################
   
  extrapolation_limits <- reactive({
    
    means <- standard_means()
    
    extrap_models <- extrapolation_models()
    
    ymin <- min(means$mean)
    ymax <- max(means$mean)
    
    # Difference between the two lowest concentrations
    low_diff <- abs(diff(extrap_models$lowest_points$mean))
    
    # Difference between the two highest concentrations
    high_diff <- abs(diff(extrap_models$highest_points$mean))
    
    list(
      ymin = ymin,
      ymax = ymax,
      
      lower_limit = ymin - input$extrapolation * low_diff,
      upper_limit = ymax + input$extrapolation * high_diff
    )
  })
  
  ###########################################
  # CALCULATE CONCENTRATIONS
  ###########################################
  calculated_data <- reactive({
    df <- raw_data()
    
    params <- model_parameters()
    
    limits <- extrapolation_limits()
    
    extrap_models <- extrapolation_models()
    
    # Empty vector
    result <- df %>% 
      mutate(calculated = NA_real_,
             model = "Out of range")
    
    # Standard curve values
    standard_curve_rows <-
      result$basic_calculation >= limits$ymin &
      result$basic_calculation <= limits$ymax
    
    result$calculated[standard_curve_rows] <-
      vapply(
        result$basic_calculation[standard_curve_rows],
        
        invert_4pl,
        
        numeric(1),
        
        lower = params$lower,
        upper = params$upper,
        ic50 = params$ic50,
        slope = params$slope
      )
    
    result$model[standard_curve_rows] <- "Standard curve"
    
    # High extrapolation
    
    high_rows <- 
      result$basic_calculation > limits$ymax &
      result$basic_calculation <= limits$upper_limit
    
    high_coef <- coef(extrap_models$highest_concentration_model)
    
    if(any(high_rows)) {
      
      result$calculated[high_rows] <- 
        (result$basic_calculation[high_rows] -
            high_coef[1]) / high_coef[2]
      
    result$model[high_rows] <-
      "Linearly extrapolated"
    }
    
    
    # Low extrapolation
    low_rows <-
      result$basic_calculation < limits$ymin &
      result$basic_calculation >= limits$lower_limit
    
    low_coef <- coef(extrap_models$lowest_concentration_model)
    
    if(any(low_rows)) {
      result$calculated[low_rows] <- 
        (result$basic_calculation[low_rows] -
            low_coef[1]) / low_coef[2]
      
      result$model[low_rows] <-
        "Linearly extrapolated"
    }
    
    result <- result %>% 
      mutate(calculated = if_else(calculated > 0,
                                  calculated, NA_real_),
             model = if_else(is.na(calculated), 
                             "Out of range", model)
             )
    
    result
    })
  
  ###########################################
  # CALCULATED CONCENTRATION LIMITS FOR PLOT
  ###########################################
  
  calc_limits <- reactive({
    calc <- calculated_data() %>% 
      filter(!is.na(calculated),
             calculated > 0)
    
    req(nrow(calc) > 0)
    
    list(xmin = min(calc$calculated),
         xmax = max(calc$calculated))
  })
  
  ###########################################
  # CREATE EXTRAPOLATION LINES FOR PLOT
  ###########################################

  extrapolation_lines <- reactive({
    limits <- calc_limits()
    
    models <- extrapolation_models()
    
    low_points <- models$lowest_points
    high_points <- models$highest_points
    
    # Low concentration line
    low_line <- NULL
    
    if(limits$xmin < min(low_points$concentration)) {
      
      low_x <- 10^seq(
        log10(limits$xmin),
        log10(min(low_points$concentration)),
        length.out = 100
      )
      
      low_line <- tibble(
        concentration = low_x,
        mean = predict(
          models$lowest_concentration_model, 
          newdata = tibble(
            concentration = low_x
          )
        )
      )
    }
    
    # High concentration line
    high_line <- NULL
    
    if(limits$xmax > max(high_points$concentration)) {
      high_x <- 10^seq(
        log10(max(high_points$concentration)),
        log10(limits$xmax),
        length.out = 100
      )
      
      high_line <- tibble(
        concentration = high_x,
        mean = predict(
          models$highest_concentration_model,
          newdata = tibble(
            concentration = high_x
          )
        )
      )
    }
    
    list(low = low_line,
         high = high_line)
  })
  
  
  ###########################################
  # OUTPUT: STANDARD TABLE
  ###########################################
  
  output$table_standard <- renderTable({
    standard_data()
    
  })
  
  ###########################################
  # OUTPUT: STANDARD CURVE
  ###########################################
  
  output$standard_curve <- renderPlot({
    fit <- model()
    calculated <- calculated_data()
    lines <- extrapolation_lines()
    
    plot <- fit$data %>% 
      ggplot(aes(x = Dose, y = Response)) +
        geom_point(size = 4, shape = 21, fill = "dodgerblue")
    
    
    if(!is.null(lines$low)) {
      
      plot <- plot +
        geom_line(
          data = lines$low,
          aes(x = concentration, y = mean), inherit.aes = FALSE
        )
    }
    
    if(!is.null(lines$high)) {
      
      plot <- plot +
        geom_line(
          data = lines$high,
          aes(x = concentration, y = mean), inherit.aes = FALSE
        )
    }
    
    plot +
      geom_point(
        data = calculated %>% filter(model %in% c("Linearly extrapolated",
                                                  "Standard curve")),
        aes(x = calculated, y = basic_calculation, fill = "model"),
        size = 3, shape = 21) +
      scale_x_log10() +
      labs(x = "Concentration", y = "Absorbance", fill = NULL) + 
      theme_classic()+
      theme(legend.position = "bottom",
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 12)
      )
  })
  
  ###########################################
  # OUTPUT: MODEL PARAMETERS
  ###########################################
  
  output$model <- renderPrint({
    model()$parameters
  })
  
  ###########################################
  # OUTPUT: CALCULATED TABLE
  ###########################################
  
  output$table_calculated  <- renderTable({
    calculated_data()
  })
  
  ###########################################
  # DOWNLOAD CSV
  ###########################################
  
  output$export_csv <- downloadHandler(
    filename = function() {
      paste0("ELISA_extrapolated_values_",
             Sys.Date(),
             ".csv")
    },
    
    content = function(file) {
      write_csv(
        calculated_data(),
        file
      )
    }
    
  )
  
  # added diagnostic inputs
  output$table_raw <- renderTable({
    raw_data()
  })
  
  output$table_metadata <- renderTable({
    concentration_data()
  })
  
  output$table_standard <- renderTable({
    standard_data()
  })
  
  
  ###########################################
  # RUN APP
  ###########################################
}


shinyApp(ui, server)
