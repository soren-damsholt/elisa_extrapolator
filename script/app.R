library(tidyverse)
library(readxl)
library(dr4pl)
library(shiny)


ui <- fluidPage(
  titlePanel("Elisa data plot and extrapolation"),
  sidebarLayout(
        # File load
        sidebarPanel(
          fileInput("file1", "Choose xlsx file", accept = ".xlsx"),
          
        # Cell ranges
          textInput(inputId = "measurement", 
                    placeholder = "e.g. A1:B2",
                    label = "Cell range of measurements"),
        
          textInput(inputId = "metadata", 
                    placeholder = "e.g. A1:B2",
                    label = "Cell range of sample data"),
        
        #Sheet name:
          textInput(inputId = "sheetname", 
                    placeholder = "e.g. Result summary",
                    label = "Sheetname with elisa summary data"),
        
        # Slider to set extrapolation factor
          sliderInput("extrapolation", label = "Extrapolation factor",
                      min = 0, max = 0.5, value = 0.1) 
          ),
        
        # Tab to go to calculated values
        
        # Button to export csv with calculated values.
    
        mainPanel(tableOutput("table_standard"),
                  plotOutput("standard_curve"),
                  verbatimTextOutput("model"),
                  tableOutput("table_calculated"))
  ),
)

  
server <- function(input, output) {
  #Helper function/calculating function
  get_concentration <- function(standard, y, lower, upper, ic50, slope, extrapolation){
    
    ymin <- min(standard$basic_calculation)
    ymax <- max(standard$basic_calculation)
    
    mean_abs <- standard %>% group_by(concentration) %>%
      summarise(mean = mean(basic_calculation), .groups = "drop")
    
    min_factor <- mean_abs %>%
      slice_min(concentration, n = 2) %>% summarise(mean_diff = mean[[2]]- mean[[1]]) %>% 
      pull(mean_diff)
    
    max_factor <- standard %>% group_by(concentration) %>%
      summarise(mean = mean(basic_calculation), .groups = "drop") %>%
      slice_max(concentration, n = 2) %>% summarise(mean_diff = mean[[1]]- mean[[2]]) %>% 
      pull(mean_diff)
    
    upper_limit <- ymax + (extrapolation * max_factor)
    lower_limit <- ymin - (extrapolation * min_factor)
    
    if(ymin <= y && y <= ymax){
      
      x = ic50 * ((((lower - upper) / (y - upper)) - 1)^(1/slope)) 
      return(x)
    }
    
    if(y > ymax && y <= upper_limit){
      
      pts <- mean_abs %>% arrange(desc(concentration)) %>% 
        slice_head(n = 2)
      
      fit <- lm(mean ~ concentration, data = pts)
      coef <- coef(fit)
      
      return((y- coef[1])/ coef[2])
    }
    
    if(y < ymin && y >= lower_limit){
      
      pts <- mean_abs %>% arrange(desc(concentration)) %>% 
        slice_tail(n = 2)
      
      fit <- lm(mean ~ concentration, data = pts)
      
      coef <- coef(fit)
      
      return((y- coef[1])/ coef[2])
      
    }
    return(NA_real_)}
     
  #Background calculations
      # Read in data and filter down to standard data
        standard_data <- reactive({
          req(input$file1,
              input$measurement, 
              input$metadata,
              input$sheetname)
          
          file <- input$file1$datapath
          
          data <- read_xlsx(file,
                            range = input$measurement,
                            sheet = input$sheetname,
                            col_names = c("well", "group", "sample",
                                          "absorbance_450", "absorbance650",
                                          "basic_calculation"))
          
          concentration <- read_xlsx(file,
                                     range = input$metadata,
                                     sheet = input$sheetname,
                                     col_names = c("sample", "concentration"))
          
          data %>% filter(str_detect(sample, "Std")) %>% 
            left_join(concentration, by = "sample") %>% 
            drop_na() %>% mutate(concentration = as.numeric(concentration))
          })
        # DR4PL model from standarad data
        model <- reactive({
          df <- standard_data()
          model <- dr4pl(
            basic_calculation ~ concentration,
            data= df,
            method.init = "logistic")
        })
        
        # read in raw values again:
        raw_data <- reactive({
          req(input$file1,
              input$measurement,
              input$sheetname)
          
          
          read_xlsx(input$file1$datapath,
                    range = input$measurement,
                    sheet = input$sheetname,
                    col_names = c("well", "group", "sample",
                                  "absorbance_450", "absorbance650",
                                  "basic_calculation"))
        })
        # Calculated values
        calculated_data <- reactive({
      df <- raw_data()
      standard <- standard_data()
      fit <- model()
      lower <- fit$parameters[4]
      upper <- fit$parameters[1]
      ic50 <- fit$parameters[2]
      slope <- fit$parameters[3]
      
      
     
      
        
    result <- df %>% rowwise() %>% 
      mutate(calculated = get_concentration(standard, basic_calculation,
                                            lower = lower, upper = upper,
                                            ic50 = ic50, slope = slope, extrapolation = input$extrapolation)) %>% 
          ungroup() %>% mutate(model = case_when(basic_calculation >= min(standard$basic_calculation) &
                                                   basic_calculation <= max(standard$basic_calculation) ~ "standard curve",
                                                 is.na(calculated) ~ "out of range",
                                                 TRUE ~ "linearly extrapolated"))
    return((result))
    
    })
        
        #
        calc_limits <- reactive({
          calc <- calculated_data() %>%
            filter(!is.na(calculated))
          req(nrow(calc) > 0)
          list(
            xmin = min(calc$calculated),
            xmax = max(calc$calculated))})
        
        # Linear rendered points for the plot
        get_linear_line_log10space <- function(standard, xmin, xmax, direction){
          
          mean_abs <- standard %>% group_by(concentration) %>% 
            summarise(mean = mean(basic_calculation), .groups = "drop")
          
        if(direction == "low"){
          pts <- mean_abs %>% arrange(desc(concentration)) %>% 
            slice_tail(n = 2)
          
          fit <- lm(mean ~ concentration, data = pts)
        
          linedata <- tibble(concentration = 10^seq(
            log10(xmin),
            log10(max(pts$concentration)),
            length.out = 100))
          
          linedata$mean <- predict(fit, newdata = linedata)
          
          return(linedata)
          
          } else{
          pts <- mean_abs %>% arrange(desc(concentration)) %>% 
            slice_head(n = 2)
          
          fit <- lm(mean ~ concentration, data = pts)
          
          linedata <- tibble(concentration = 10^seq(
            log10(min(pts$concentration)),
            log10(xmax),
            length.out = 100))
          
          linedata$mean <- predict(fit, newdata = linedata)
          
          return(linedata)
          }
        }
        
      lowline <- reactive({
        limits <- calc_limits()
        
        get_linear_line_log10space(standard = standard_data(),
                                    xmin = limits$xmin,
                                    xmax = limits$xmax,
                                    direction = "low")})
      
      highline <- reactive({
        limits <- calc_limits()
        
        get_linear_line_log10space(standard = standard_data(),
                                    xmin = limits$xmin,
                                    xmax = limits$xmax,
                                    direction = "high")})
      
  #Rendered output
        
  output$table_standard <- renderTable({
    standard_data()
  })
    
    output$standard_curve <- renderPlot({
    fit <- model()
    lower <- fit$parameters[4]
    upper <- fit$parameters[1]
    ic50 <- fit$parameters[2]
    slope <- fit$parameters[3]
    
    calculated <- calculated_data()
    
    
    fit$data %>% ggplot(aes(x = Dose, y = Response)) + 
      geom_point(size = 5, shape = 21, fill = "dodgerblue") + 
      geom_smooth(method = "nls", se = FALSE, span = 0.2,
                  formula = y ~ lower + (upper-lower)/(1+(ic50/(x^10))^slope),
                  method.args = list(start = list(lower = lower, 
                                                  upper = upper,
                                                  ic50 = ic50,
                                                  slope = slope)),
                  fullrange = TRUE)+
      geom_line(data = lowline(), aes(concentration, mean), inherit.aes = FALSE)+
      geom_line(data = highline(), aes(concentration, mean), inherit.aes = FALSE)+
      geom_point(data = calculated %>% filter(model %in% c("linearly extrapolated",
                                                                "standard curve")), 
                 aes(x = calculated, y = basic_calculation, fill = model), size = 3, shape = 21)+
      scale_x_log10(guide = guide_axis_logticks(short = 1, mid = 1.5, long = 2))+
      coord_cartesian()+
      labs(x = "Concentration", y = "Absorbance", fill = NULL)+
      theme_classic()+
      theme(
        legend.position = "bottom",
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(5, "pt"))
    })
  
  output$model <- renderPrint({
    
    df <- model()
    
    df$parameters
    
  })
  
  output$table_calculated <- renderTable({
    calculated_data()
  })
}


shinyApp(ui, server)
