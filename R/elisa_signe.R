library(tidyverse)
library(readxl)
library(dr4pl)

data <- read_xlsx("elisa_data.xlsx", range = c("B6:G101"), sheet = "Result summary", 
                  col_names = c("well", "group", "sample", "absorbance_450", "absorbance650",
                                "basic_calculation")) 

concentration <- read_xlsx("elisa_data.xlsx", range = c("B140:C144"), sheet = "Result summary",
                           col_names = c("sample", "concentration")
                           )

standard_data <- data %>% filter(str_detect(sample, "Std")) %>% left_join(., concentration, by = "sample") %>% 
  drop_na %>% mutate(concentration = as.numeric(concentration))


model <- standard_data %>% dr4pl(basic_calculation ~ concentration, method.init = "logistic", data=.)

model %>% plot()

model_parameters <- model$parameters

# x = c * ((a - d / y - d) - 1)^1/b 
# a = theoretical response at 0 concentration. b = slope factor,
# c = mid-range concentration (inflection point), d = theoretical concentration at infinite concentration

get_concentration <- function(standard, y, parameters, extrapolation = 0.10){
  
  ymin <- min(standard$basic_calculation)
  ymax <- max(standard$basic_calculation)
  
  mean_abs <- standard %>% group_by(concentration) %>%
    summarise(mean = mean(basic_calculation), .groups = "drop")
  
  min_factor <- mean_abs %>%
    slice_min(concentration, n = 2) %>% summarise(mean_diff = mean[[2]]- mean[[1]])
  
  max_factor <- standard %>% group_by(concentration) %>%
    summarise(mean = mean(basic_calculation), .groups = "drop") %>%
    slice_max(concentration, n = 2) %>% summarise(mean_diff = mean[[1]]- mean[[2]])
  
  upper_limit <- ymax + (extrapolation * max_factor)
  lower_limit <- ymin - (extrapolation * min_factor)
  
  if(ymin <= y && y <= ymax){
  
  x = parameters[2] * ((((parameters[4] - parameters[1]) / (y - parameters[1])) - 1)^(1/parameters[3])) 
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
  return(NA_real_)
  
}


calculated_data <- data %>% rowwise() %>% mutate(calculated = get_concentration(standard_data, basic_calculation,
                                                            model_parameters, extrapolation = 0.5)) %>% 
  ungroup() %>% mutate(model = case_when(basic_calculation >= min(standard_data$basic_calculation) &
                                           basic_calculation <= max(standard_data$basic_calculation) ~ "standard curve",
                                         is.na(calculated) ~ "out of range",
                                         TRUE ~ "linearly extrapolated"))



get_linear_line_log10space <-function(standard, direction){
    
    if(direction == "low"){
      
    lowline <- standard_data %>% group_by(concentration) %>% 
      summarise(mean = mean(basic_calculation), .groups = "drop") %>% 
      arrange(desc(concentration)) %>% slice_tail(n = 2) 
    
    lowline_model <- lowline %>% 
      lm(mean~concentration, data=.)
    
    linedata <- tibble(concentration = 10^seq(
                      log10(min(calculated_data$calculated, na.rm = TRUE)),
                      log10(max(lowline$concentration)),
                      length.out = 100)) %>% 
      mutate(
        mean = predict(lowline_model, 
        newdata = data.frame(concentration = concentration)))
    
return(linedata)
    }
  
  if(direction == "high"){
    
    highline <- standard %>% group_by(concentration) %>% 
      summarise(mean = mean(basic_calculation), .groups = "drop") %>% 
      arrange(desc(concentration)) %>% slice_head(n = 2) 
  
    highline_model <- highline %>% 
      lm(mean~concentration, data=.)
  
    linedata <- tibble(concentration = 10^seq(
      log10(min(highline$concentration)),
      log10(max(calculated_data$calculated, na.rm = TRUE)),
      length.out = 100)) %>% 
      mutate(
        mean = predict(highline_model, 
                       newdata = data.frame(concentration = concentration)))
    
  return(linedata)}
    
}

lowline <- get_linear_line_log10space(standard_data, "low")
highline <- get_linear_line_log10space(standard_data, "high")

# Manually plotting the curve in log10 space
lower <- model_parameters[4]
upper <- model_parameters[1]
ic50 <- model_parameters[2]
slope <- model_parameters[3]

model$data %>% ggplot(aes(x = Dose, y = Response)) + 
  geom_point(size = 5, shape = 21, fill = "dodgerblue") + 
  geom_smooth(method = "nls", se = FALSE, span = 0.2,
              formula = y ~ lower + (upper-lower)/(1+(ic50/(x^10))^slope),
              method.args = list(start = list(lower = lower, 
                                              upper = upper,
                                              ic50 = ic50,
                                              slope = slope)),
              fullrange = TRUE)+
  geom_line(data = lowline, aes(concentration, mean))+
  geom_line(data = highline, aes(concentration, mean))+
  geom_point(data = calculated_data %>% filter(model %in% c("linearly extrapolated",
                                                            "standard curve")), 
             aes(x = calculated, y = basic_calculation, fill = model), size = 3, shape = 21)+
  scale_x_log10(guide = guide_axis_logticks(short = 1, mid = 1.5, long = 2))+
  coord_cartesian()+
  labs(x = "Concentration", y = "Absorbance", fill = NULL)+
  theme_classic()+
  theme(
    legend.position = "bottom",
  )
