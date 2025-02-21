library(imfr)
databases <- imfr::imf_databases('IFS') %>% 
  as_tibble()

imfr::imf_codelist('IFS')
x <- imfr::imf_parameters('IFS')[['indicator']] %>% 
  as_tibble() %>% 
  filter(str_detect(description %>% tolower, 'rate'))

x %>% 
  filter(str_detect(description %>% tolower, 'year'))
  print(n = 200)
"
90 FIMM_PA                       Financial, Interest Rates, Money Market, Percent per annum                           
86 FPOLM_PA                      Financial, Interest Rates, Monetary Policy-Related Interest Rate, Percent per annum  
82 FILIBOR_1M_PA                 Financial, Interest Rates, London Interbank Offer Rate, One-Month, Percent per annum 
49 FITB_3M_PA                    Financial, Interest Rates, Government Securities, Treasury Bills, 3-month, Percent p~
30 FII_3M_PA                     Financial, Interest Rates, 3-Month Interbank Interest, Percent per annum             
"

ir <- imfr::imf_dataset('ifs', 
                  indicator = c('FIMM_PA', 'FPOLM_PA', 'FILIBOR_1M_PA', 'FITB_3M_PA', 'FII_3M_PA'),
                  freq = 'Q')

ir <- ir %>% 
  as_tibble() %>% 
  select(ccode2 = ref_area, indicator, date, value) %>% 
  mutate(year = substr(date, 1, 4),
         month = case_when(str_detect(date, 'Q1') ~ '01',
                           str_detect(date, 'Q2') ~ '04',
                           str_detect(date, 'Q3') ~ '07',
                           str_detect(date, 'Q4') ~ '10'),
         date = paste(year, month, '01', sep = '-') %>% as_date(),
         value = as.numeric(value)) %>% 
  select(-year, -month) %>% 
  group_by(ccode2, indicator) %>% 
  drop_na() %>% 
  mutate(value_fd = value - lag(value))




ir <- ir %>% 
  nest() %>%
  mutate(data = map(data, ~tryCatch({.x %>% 
                      mutate(value_ar = value %>% 
                               ts() %>% 
                               arima(order = c(4,0,0)) %>% 
                               resid() )},
                      error = function(error){NA} ))) %>% 
  unnest(data)

ir <- ir %>% 
  #Create RHS
  mutate(map_dfc(seq(4), ~ lag(value_fd, n = .x)) %>%
           set_names(paste('shockFD_lag', seq(4),sep = ''))) %>% 
  #Create RHS
  mutate(map_dfc(seq(4), ~ lag(value_ar, n = .x)) %>%
           set_names(paste('shockAR_lag', seq(4),sep = ''))) %>% 
  rename(shockFD = value_fd,
         shockAR = value_ar) %>% 
  nest()
  
ir <- ir %>% 
  unnest(data) %>% 
  ungroup(ccode2) %>% 
  nest()

ir <- ir %>% 
  unnest(data) %>% 
  group_by(ccode2, indicator) %>% 
  mutate(n = n()) %>% 
  ungroup(indicator) %>% 
  filter(n == max(n)) 
