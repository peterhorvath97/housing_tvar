data.query <- function(){
  residential <- function(){
    require(tidyverse)
    require(BISdata)
    require(lubridate)
    require(stringr)
    require(rebus)
    
    data <- BISdata::fetch_dataset(dest.dir = tempdir(),
                                   'WS_SPP_csv_col.zip') %>% 
      as_tibble %>% 
      gather(key = date, 
             value = value, 
             contains('\\d\\d\\d\\d-Q\\d'))
    
    
    data <- data %>% 
      gather(key = date, 
             value = value, 
             contains('-Q')) %>% 
      filter(VALUE == 'R',
             `Unit of measure` == 'Index, 2010 = 100') %>% 
      select( 
        ccode2 = REF_AREA,
        date, value) %>% 
      arrange(ccode2, date) %>% 
      drop_na() %>%
      mutate(year = substr(date, 1, 4),
             month = case_when(str_detect(date, 'Q1') ~ '01',
                               str_detect(date, 'Q2') ~ '04',
                               str_detect(date, 'Q3') ~ '07',
                               str_detect(date, 'Q4') ~ '10'),
             date = paste(year, month, '01', sep = '-') %>% as_date()) %>% 
      select(-year, -month) %>% 
      rename(RPPI = value)
    
    data
    
  }
  
  intrates <- function(){
    require(tidyverse)
    require(imfr)
    require(lubridate)
    require(stringr)
    
    data <- imfr::imf_dataset('ifs', 
                              #ref_area = countries,
                              indicator = c('FIMM_PA', 'FPOLM_PA', 'FILIBOR_1M_PA', 'FITB_3M_PA', 'FII_3M_PA'),
                              freq = 'Q')
    
    data <- data %>% 
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
      rename(R = value)
    
    
    data <- data %>% 
      group_by(ccode2, indicator) %>% 
      mutate(n = n()) %>% 
      ungroup(indicator) %>% 
      filter(n == max(n)) %>% 
      ungroup() %>% 
      select(-n)
    
    
    data <- data %>% 
      select(ccode2, date, R)
    
    
    
    data2 <- bind_rows(fredr('IR3TIB01ATM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'AT'),
                       fredr('IR3TIB01BEM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'BE'),
                       fredr('IR3TIB01GRM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'GR'),
                       fredr('LVAIR3TIB01STM', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'LV'),
                       fredr('IR3TIB01SKM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'SK'),
                       fredr('IR3TIB01LUM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'LU'),
                       fredr('IR3TIB01EEM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'EE'),
                       fredr('IR3TIB01PTM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'PT'),
                       fredr('LTUIR3TIB01STM', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'LT'),
                       fredr('IRSTCI01CHM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'CH'),
                       fredr('IR3TIB01NLM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'NL'),
                       fredr('IR3TIB01FRM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'FR'),
                       fredr('IR3TIB01DEM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'DE'),
                       fredr('IR3TIB01GBM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'GB'),
                       fredr('IR3TIB01IDM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'ID'),
                       fredr('IR3TIB01SIM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'SI'),
                       fredr('IRSTCI01JPM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'JP'),
                       fredr('IR3TIB01IEM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'IE'),
                       fredr('IR3TIB01CAM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'CA'),
                       fredr('IR3TIB01SEM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'SE'),
                       fredr('IR3TIB01FIM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'FI'),
                       fredr('IR3TIB01ITM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'IT'),
                       fredr('IR3TIB01ESM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'ES'),
                       fredr('FEDFUNDS', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'US'),
                       fredr('IRSTCI01INM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'IN'),
                       fredr('IR3TIB01IDM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'ID'),
                       fredr('IR3TIB01PLM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'PL'),
                       fredr('INTDSRTRM193N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'TR'),
                       fredr('IR3TIB01CZM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'CZ'),
                       fredr('COLIR3TIB01STM', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'CO'),
                       fredr('IRSTCI01MXM156N', frequency = 'q') %>% 
                         drop_na() %>% 
                         select(date, R = value) %>% 
                         mutate(ccode2 = 'MX')
                       
    )
    
    
    data <- data %>% 
      filter(!(ccode2 %in% data2$ccode2)) %>% 
      bind_rows(data2)
    
    data
    
  }
  
  gdp <- function(){
    require(tidyverse)
    require(imfr)
    require(lubridate)
    require(stringr)
    require(fredr)
    
    
    data1 <- imfr::imf_dataset('ifs', 
                               #ref_area = countries[i],
                               indicator = 'NGDP_SA_XDC',
                               freq = 'Q')
    
    data2 <- imfr::imf_dataset('ifs', 
                               #ref_area = countries[i],
                               indicator = 'NGDP_D_SA_IX',
                               freq = 'Q')
    
    
    data1 <- data1 %>% 
      as_tibble() %>% 
      select(ccode2 = ref_area, date, value) %>% 
      mutate(year = substr(date, 1, 4),
             month = case_when(str_detect(date, 'Q1') ~ '01',
                               str_detect(date, 'Q2') ~ '04',
                               str_detect(date, 'Q3') ~ '07',
                               str_detect(date, 'Q4') ~ '10'),
             date = paste(year, month, '01', sep = '-') %>% as_date(),
             value = as.numeric(value)) %>% 
      select(-year, -month) %>% 
      rename(Y = value)
    
    
    data2 <- data2 %>% 
      as_tibble() %>% 
      select(ccode2 = ref_area, date, value) %>% 
      mutate(year = substr(date, 1, 4),
             month = case_when(str_detect(date, 'Q1') ~ '01',
                               str_detect(date, 'Q2') ~ '04',
                               str_detect(date, 'Q3') ~ '07',
                               str_detect(date, 'Q4') ~ '10'),
             date = paste(year, month, '01', sep = '-') %>% as_date(),
             value = as.numeric(value)) %>% 
      select(-year, -month) %>% 
      rename(P = value)
    
    data <- data1 %>% 
      inner_join(data2) %>% 
      mutate(Y = 100*Y/P) %>% 
      select(-P)
    
    
    
    data <- data %>% 
      bind_rows(fredr('NAEXKP01RUQ652S') %>% 
                  select(date, Y = value) %>% 
                  mutate(ccode2 = 'RU')) %>% 
      bind_rows(fredr('CLVMNACSAB1GQIS') %>% 
                  select(date, Y = value) %>% 
                  mutate(ccode2 = 'IS'))
    
    
    
    data <- data %>% 
      mutate(year = year(date)) %>% 
      group_by(ccode2, year) %>% 
      mutate(mean = mean(Y),
             Y2010 = ifelse(year == 2010,1,NA),
             mean = mean*Y2010) %>% 
      ungroup(year) %>% 
      mutate(mean = mean(mean, na.rm = T),
             Y = 100*Y/mean) %>% 
      ungroup()
    
    
    data <- data %>% 
      select(ccode2, date, Y)
    
    data
  }
  
  cpi <- function(){
    require(tidyverse)
    require(imfr)
    require(lubridate)
    require(stringr)
    
    data <-   imfr::imf_dataset('ifs', 
                                #ref_area = countries,
                                indicator = 'PCPI_IX',
                                freq = 'Q')
    
    
    
    data <- data %>% 
      as_tibble() %>% 
      select(ccode2 = ref_area, date, value) %>% 
      mutate(year = substr(date, 1, 4),
             month = case_when(str_detect(date, 'Q1') ~ '01',
                               str_detect(date, 'Q2') ~ '04',
                               str_detect(date, 'Q3') ~ '07',
                               str_detect(date, 'Q4') ~ '10'),
             date = paste(year, month, '01', sep = '-') %>% as_date(),
             value = as.numeric(value)) %>% 
      select(-year, -month) %>% 
      rename(P = value)
    
    data
    
  }
  
  
  data_residential <- residential()
  #countries <- unique(data_residential$ccode2)
  
  data_intrates <- intrates()
  #Sys.sleep(runif(1,2,6))
  #countries <- unique(data_intrates$ccode2)
  #countries <- countries[!(countries %in% c('MA', 'PE', 'CN', 'IS', 'RU', 'MY'))]
  
  data_gdp <- gdp()
  #Sys.sleep(runif(1,2,6))
  data_cpi <- cpi()
  #Sys.sleep(runif(1,2,6))
  #Sys.sleep(runif(1,2,6))
  
  
  
  dbnames <- ls()[str_detect(ls(), 'data_')]
  
  data <- get(dbnames[1])
  
  for(name in dbnames[-1]){
    data <- full_join(data, get(name))
  }
  
  
  
  data <- data %>% 
    filter(!is.na(P),
           !is.na(Y),
           !is.na(R))
  
  data
}

data <- data.query()