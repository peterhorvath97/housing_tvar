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


data <- fred_query(ids = 
                 c('QUSR628BIS',
                   'DFF',
                   'MORTGAGE30US',
                   'PCECC96',
                   'REALLN',
                   'FEDFUNDS',
                   'GDPC1',
                   'PCEPI',
                   'GS10', 'GS1',
                   'BAA', 'AAA'
                   ),
           'q')


data <- data %>% 
  mutate(R = DFF,
         MR = MORTGAGE30US,
         CS = BAA-AAA,
         TS = GS10 - GS1,
         HL = REALLN,
         HP = QUSR628BIS,
         P = PCEPI,
         Y = GDPC1) %>% 
  select(date, R, MR, CS, TS, HL, HP, P, Y) %>% 
  mutate(Y = log(Y), 
         P = log(P),
         HP = log(HP),
         HL = log(HL))


shadow <- readxl::read_excel('data/shadowrate_US.xls', col_names = F) %>% 
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


data <- data %>% 
  full_join(shadow, by = 'date') %>% 
  mutate(R = ifelse(round(R.x) == 0, R.y, R.x)) %>% 
  select(-R.x, -R.y)
rm(shadow)


data <- rdbnomics::rdb('BIS',
               'WS_CREDIT_GAP',
               'Q.US.P.A.C') %>% 
  as_tibble() %>% 
  select(date = period, CGG = value) %>% 
  inner_join(data)


rent_price <- fred_query(ids = 
                     c('QUSR628BIS',
                       'CUUR0000SEHA',
                       'CUUR0000SA0L2'),
                   'q')

rent_price <- rent_price %>% 
  transform_to_base_index(c('CUUR0000SEHA', 'QUSR628BIS', 'CUUR0000SA0L2'), 2015) %>% 
  mutate(CUUR0000SEHA = 100*CUUR0000SEHA/CUUR0000SA0L2,
         RP = log(CUUR0000SEHA) - log(QUSR628BIS)) %>% 
  select(date, RP) 



data <- inner_join(data, rent_price)


data <- data %>% 
  mutate(RP = scale(RP)[,1],
         CGG = scale(CGG)[,1])



p <- data %>% 
  gather(var, value, -date) %>% 
  mutate(var = case_when(var == 'CS' ~ 'Credit Spread',
                         var == 'HP' ~ 'Property Prices',
                         var == 'P' ~ 'Consumer Prices',
                         var == 'TS' ~ 'Term Spread',
                         var == 'HL' ~ 'Real Estate Loans',
                         var == 'MR' ~ 'Mortgage Rate',
                         var == 'R' ~ 'Fed Funds Rate',
                         var == 'Y' ~ 'Real GDP',
                         var == 'CGG' ~ 'Credit-to-GDP gap',
                         var == 'RP' ~ 'Rent-price ratio')) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_line(linewidth = 1) +
  facet_wrap(~var, scales = 'free', ncol = 2) +
  theme_minimal() +
  labs(x = '', y = '')
ggsave('figures/data_tsplots.pdf', p)


data

}
