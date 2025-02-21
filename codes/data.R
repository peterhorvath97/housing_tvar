data <- function(){
# Function to transform selected variables to a base year index
transform_to_base_index <- function(data, variables = NULL, base_year) {
  
  require(dplyr)
  
  # Ensure the input data has a date column
  if (!"date" %in% names(data)) {
    stop("The dataset must contain a 'date' column.")
  }
  
  # Set default variables to all columns except the date column
  if (is.null(variables)) {
    variables <- setdiff(names(data), "date")
  }
  
  # Extract the average values for the base year
  base_row <- data %>%
    filter(format(date, "%Y") == as.character(base_year)) %>%
    summarize(across(all_of(variables), \(x) mean(x, na.rm = TRUE)))
  
  if (nrow(base_row) == 0) {
    stop("No data found for the specified base year.")
  }
  
  # Extract the base year values for the selected variables
  base_values <- base_row %>%
    select(all_of(variables)) %>%
    unlist()
  
  # Transform the selected variables
  data_transformed <- data %>%
    mutate(across(all_of(variables), ~ . / base_values[deparse(substitute(.))] * 100))
  
  return(data_transformed)
}
fred_query <- function(ids, freq, start = '1940-01-01'){
  
  require(fredr)
  require(dplyr)
  require(purrr)
  
  #download data
  params <- list(
    series_id = ids,
    frequency = freq,
    observation_start = as.Date(start)
  )
  
  
  data  <- pmap_dfr(
    .l = params,
    .f = ~ fredr(series_id = .x, frequency = .y) ) %>%
    dplyr::select(date, series_id, value) %>%
    spread(key = series_id, value = value) %>%
    drop_na() 
  
  data
}
x13adj <- function(value, freq, date){
  require(forecast)
  require(seasonal)
  require(tidyverse)
  require(lubridate)
  tryCatch(
    value %>% 
      ts(frequency = freq,
         start = c(date %>% min %>% year, 
                   date %>% min %>% month)) %>% 
      seas() %>% 
      seasadj() %>% 
      as.numeric(),
    error = function(error){value}
  )
}

fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

# House & rent prices block
data <- fredr('CSUSHPISA') %>% 
  select(date, HP = value) %>% 
  drop_na() %>% 
  transform_to_base_index(base_year = 2015) %>% 
  inner_join(
    fredr('CUSR0000SEHA') %>% 
      select(date, R = value) %>% 
      drop_na() %>% 
      transform_to_base_index(base_year = 2015)
    
  ) %>% 
  inner_join(
    fredr('CUUR0000SA0L2') %>% 
      select(date, P = value) %>% 
      drop_na() %>% 
      transform_to_base_index(base_year = 2015)
  ) %>% 
  mutate(HP = 100*HP/P, 
         R = 100*R/P) %>% 
  select(-P) 



data %>% 
  gather(var, value, -date) %>% 
  ggplot(aes(x = date, y = value, color = var)) +
  geom_line(linewidth = 1)


data %>% 
  gather(var, value, -date) %>% 
  group_by(var) %>% 
  mutate(value = value-lag(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = date, y = value, color = var)) +
  geom_line(linewidth = 1)


data %>% 
  mutate(ret = (R+HP-lag(HP))/lag(HP)) %>% 
  ggplot(aes(x = date, y = ret)) +
  geom_line(linewidth = 1)


data %>% 
  mutate(rp = log(R) - log(HP)) %>% 
  ggplot(aes(x = date, y = rp)) +
  geom_line(linewidth = 1)

data <- data %>% 
  mutate(RP = log(R) - log(HP))


# Real econ block
macro <- fred_query(ids = c('FEDFUNDS', 'CPIAUCSL', 'INDPRO'), freq = 'm') %>% 
  rename(P = CPIAUCSL,
         i = FEDFUNDS,
         Y = INDPRO) %>% 
  transform_to_base_index(c('P', 'Y'), 2015)

# Consumption

cons <- fred_query(ids = c('PCE', 'PCEPI'), freq = 'm') %>% 
  transform_to_base_index(variables = 'PCEPI', 2015) %>% 
  mutate(C = 100*PCE/PCEPI) %>% 
  select(date, C)

# Mortgage market block

mort <- fred_query(ids = c('MORTGAGE30US', 'RHEACBW027SBOG'), freq = 'm') %>% 
  rename(MR = MORTGAGE30US,
         HL = RHEACBW027SBOG) 

# Exogenous regressors

exog <- fred_query(ids = c('POPTHM', 'PAYEMS'), freq = 'm') %>% 
  rename(POP = POPTHM,
         EMP = PAYEMS)

out <- data %>% 
  inner_join(macro) %>% 
  inner_join(mort) %>% 
  inner_join(cons) %>% 
  inner_join(exog) %>% 
  mutate(across(c(HP, R, P, Y, HL, C, EMP, POP), ~100*.x/lag(.x,12)-100 )) %>% 
  drop_na() 


out %>% 
  gather(key, value, -date) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~key, scales = 'free')

out
}
