install.packages('BISdata')
library(tidyverse)
library(lubridate)
library(BISdata)
datasets <- BISdata::datasets()

files <- datasets %>% 
  filter(str_detect(filename %>% tolower, 'cpp')) %>% 
  pull(filename)

files <- c(
'WS_SPP_csv_col.zip',
'WS_CREDIT_GAP_csv_col.zip'
)

data <- vector(mode = 'list', length = length(files))
names(data) <- files

for(i in 1:length(data)){
  data[[i]] <- fetch_dataset(dest.dir = tempdir(),
                             files[i])
}

lapply(data, as_tibble %>% 
         gather(key = date, 
                value = value, 
                contains('\\d\\d\\d\\d-Q\\d')))


hprice <- data[[1]] %>% 
  as_tibble %>% 
  gather(key = date, 
         value = value, 
         contains('-Q')) %>% 
  filter(VALUE == 'R',
         `Unit of measure` == 'Index, 2010 = 100') %>% 
  select( 
         ccode2 = REF_AREA,
         country = `Reference area`,
         date, value) %>% 
  arrange(country, date) %>% 
  drop_na() %>%
  mutate(year = substr(date, 1, 4),
         month = case_when(str_detect(date, 'Q1') ~ '01',
                           str_detect(date, 'Q2') ~ '04',
                           str_detect(date, 'Q3') ~ '07',
                           str_detect(date, 'Q4') ~ '10'),
         date = paste(year, month, '01', sep = '-') %>% as_date()) %>% 
  select(-year, -month) 

#creditgap maybe

hprice <- hprice %>% 
  group_by(country) %>% 
  #mutate(lvalue = value - lag(value, 4)) %>% 
  #mutate(lvalue = value - lag(value)) %>% 
  mutate(lvalue = value) %>% 
  drop_na() %>% 
  mutate(value = value - mean(value),
         ind = case_when(lvalue <= quantile(lvalue, .1) ~ 'q1',
                         lvalue <= quantile(lvalue, .2) ~ 'q2',
                         lvalue <= quantile(lvalue, .3) ~ 'q3',
                         lvalue <= quantile(lvalue, .4) ~ 'q4',
                         lvalue <= quantile(lvalue, .5) ~ 'q5',
                         lvalue <= quantile(lvalue, .6) ~ 'q6',
                         lvalue <= quantile(lvalue, .7) ~ 'q7',
                         lvalue <= quantile(lvalue, .8) ~ 'q8',
                         lvalue <= quantile(lvalue, .9) ~ 'q9',
                         lvalue <= quantile(lvalue, 1) ~ 'q10')) %>% 
  select(-lvalue) %>% 
  #Create RHS
  mutate(map_dfc(seq(4), ~ lag(value, n = .x)) %>%
           set_names(paste('lag', seq(4),sep = ''))) %>% 
  #Create LHS
  mutate(map_dfc(seq(21), ~ lead(value, n = .x)) %>%
           set_names(paste('lead', seq(21),sep = '')))
