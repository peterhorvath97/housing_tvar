data <- function(){
  #import data
  #set api key
  fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")
  
  
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
  
  seasonal_adjustment <- function(data, date_col) {
    # Ensure the date column is properly evaluated
    date_col <- ensym(date_col)
    
    # Convert date column to date format if not already
    data <- data %>% mutate(!!date_col := as.Date(!!date_col))
    
    # Extract variable names (exclude the date column)
    vars <- data %>% select(-!!date_col) %>% names()
    
    # Compute frequency based on the date intervals
    date_diff <- median(diff(sort(pull(data, !!date_col))))
    frequency <- case_when(
      date_diff >= 360 ~ 1,  # Annual data
      date_diff >= 85 & date_diff < 360 ~ 4,  # Quarterly data
      date_diff >= 28 & date_diff < 85 ~ 12,  # Monthly data
      TRUE ~ 52  # Weekly or more granular
    )
    
    message("Detected frequency: ", frequency)
    
    # Create a list to store seasonally adjusted variables
    adjusted_data <- data
    
    # Loop through variables and check seasonality
    for (var in vars) {
      message("Processing variable: ", var)
      
      # Create a time series object
      ts_data <- ts(data[[var]], start = c(year(min(pull(data, !!date_col))), 
                                           month(min(pull(data, !!date_col)))),
                    frequency = frequency)
      
      # Check if the series has enough data points for seasonality analysis
      min_data_points <- frequency * 2  # Require at least two full years
      if (length(ts_data) < min_data_points) {
        message("Skipping ", var, ": not enough data.")
        next
      }
      
      # Try seasonal adjustment
      tryCatch({
        seas_result <- seas(ts_data)
        
        # Check if the seasonal component is significant
        if (!is.null(seas_result$estimates$seas) && any(seas_result$estimates$seas != 0)) {
          message("Seasonality detected in ", var, ". Seasonally adjusting...")
          adjusted_series <- final(seas_result)
          
          # Replace the variable with its seasonally adjusted counterpart
          adjusted_data[[var]] <- as.numeric(adjusted_series)
        } else {
          message("No significant seasonality detected for ", var)
        }
      }, error = function(e) {
        message("Error adjusting variable ", var, ": ", e$message)
      })
    }
    
    return(adjusted_data)
  }
  
  
  read_excel('data/shadowrate_US.xls', col_names = F) %>% 
    rename(date = ...1,
           R = ...2) %>% 
    mutate(date = paste(
      substr(date, 1,4), 
      substr(date, 5,6),
      '01', sep = '-') %>% as_date()) %>% 
    mutate(year = year(date),
           quarter = quarter(date)) %>% 
    group_by(year, quarter) %>% 
    mutate(R = mean(R)) %>% 
    ungroup() %>% 
    mutate(date = paste(
      year(date),
      case_when(quarter == 1 ~ '01',
             quarter == 2 ~ '04',
             quarter == 3 ~ '07',
             quarter == 4 ~ '10'),
      '01', sep = '-') %>% as_date()) %>% 
    distinct(date, R)
      
      
      dataQ <- fred_query(
        ids = c('QUSR628BIS', #HPRICE
                'RHEACBW027SBOG', #HLOAN
                'MORTGAGE30US', #M30
                'DFF', #R
                'CPIAUCSL', #P
                'GDPC1' #Y
        ),
        freq = 'q'
      ) %>% 
        rename(P = CPIAUCSL,
               R = DFF,
               Y = GDPC1,
               M30 = MORTGAGE30US,
               HPRICE = QUSR628BIS,
               HLOAN = RHEACBW027SBOG) %>% 
        select(date, R, M30, HLOAN, HPRICE, P, Y) %>% 
        transform_to_base_index(variables = c('P', 'Y', 'HPRICE', 'HLOAN'),
                                base_year = 2015) %>% 
        mutate(P = P - lag(P), 
               Y = Y - lag(Y),
               HPRICE = HPRICE - lag(HPRICE),
               HLOAN = HLOAN - lag(HLOAN)) %>% 
        drop_na()
      
      dataQplus  <- fred_query(
        ids = c('QUSR628BIS', #HPRICE
                'RHEACBW027SBOG', #HLOAN
                'MORTGAGE30US', #M30
                'DFF', #R
                'CPIAUCSL', #P
                'GDPC1', #Y
                'PCEC', #C
                'GPDIC1' #Inv
        ),
        freq = 'q'
      ) %>% 
        rename(P = CPIAUCSL,
               R = DFF,
               Y = GDPC1,
               M30 = MORTGAGE30US,
               HPRICE = QUSR628BIS,
               HLOAN = RHEACBW027SBOG,
               C = PCEC,
               Inv = GPDIC1) %>% 
        select(date, R, M30, HLOAN, HPRICE, P, C, Inv, Y) %>% 
        transform_to_base_index(variables = c('P', 'Y', 'HPRICE', 'HLOAN', 'C', 'Inv'),
                                base_year = 2015) %>% 
        mutate(P = P - lag(P), 
               Y = Y - lag(Y),
               HPRICE = HPRICE - lag(HPRICE),
               HLOAN = HLOAN - lag(HLOAN),
               C = C - lag(C),
               Inv = Inv - lag(Inv)) %>% 
        drop_na()
      
      
      
      dataM <- fred_query(
        ids = c('CSUSHPISA', #HPRICE
                'RHEACBW027SBOG', #HLOAN
                'MORTGAGE30US', #M30
                'DFF', #R
                'CPIAUCSL', #P
                'INDPRO' #IP
        ),
        freq = 'm'
      )  %>% 
        rename(P = CPIAUCSL,
               R = DFF,
               Y = INDPRO,
               M30 = MORTGAGE30US,
               HPRICE = CSUSHPISA,
               HLOAN = RHEACBW027SBOG) %>% 
        select(date, R, M30, HLOAN, HPRICE, P, Y) %>% 
        transform_to_base_index(variables = c('P', 'Y', 'HPRICE', 'HLOAN'),
                                base_year = 2015) %>% 
        mutate(HPRICE = 100*HPRICE/P) %>% 
        mutate(P = P - lag(P), 
               Y = Y - lag(Y),
               HPRICE = HPRICE - lag(HPRICE),
               HLOAN = HLOAN - lag(HLOAN)) %>% 
        drop_na()
      
      
      out <- list(dataQ = dataQ,
                  dataQplus = dataQplus,
                  dataM = dataM)
      
      
}