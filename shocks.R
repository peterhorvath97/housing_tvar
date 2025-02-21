est_shocks <- function(data){
paste('codes/varsignr/', list.files('codes/varsignr') , sep = '') %>% lapply(source)


shocks <- data %>% 
  select(-RPPI) %>% 
  select(-date) %>% 
  group_by(ccode2) %>% 
  nest() %>% 
  #filter(ccode2 == 'US') %>% 
  mutate(data = map(data, ~.x %>% 
                      select(R, P, Y)),
         data = map(data, ~.x %>% 
                      ts()),
         model = map(data, ~.x %>% rwz.reject(constrained = c(+1,-2,-3), KMIN = 4, KMAX = 8)),
         shocks = map(model, ~.x[['SHOCKS']] %>% 
                        as_tibble() %>% 
                        gather(key, value, everything()) %>% 
                        group_by(key) %>% 
                        mutate(value = median(value)) %>% 
                        ungroup() %>% 
                        distinct(key, .keep_all = T) %>% 
                        mutate(key = row_number()))) %>% 
  select(-model, -data)
  
shocks <- shocks %>% 
  unnest(shocks) %>% 
  select(-key) %>% 
  rename(shock = value)




data <- data %>% 
  group_by(ccode2) %>% 
  nest() %>% 
  mutate(data = map(data, ~.x %>% 
                      slice(5:nrow(.x)))) %>% 
  unnest(data) %>% 
  ungroup() %>% 
  bind_cols(shocks %>% ungroup %>%  select(shock)) %>% 
  drop_na()


data
}

data_shocks <- est_shocks(data)

