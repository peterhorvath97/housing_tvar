#Split data
data_train <- data[1:(nrow(data)-12), ]
data_test <- data %>% 
  anti_join(data_train, by = 'date')


# Models
models <- NULL

#AR1
models[['AR1']] <- data_train %>% 
  select(HP) %>% 
  ts() %>% 
  Arima(order = c(1, 0, 0))


#ARiMA optimal
models[['ARiMA']] <- data_train %>% 
  select(HP) %>% 
  ts() %>% 
  auto.arima()


#VAR
models[['VAR']] <- data_train %>% 
  select(-date) %>% 
  ts() %>% 
  vars::VAR(p = 1, type = 'const')


#TVAR1
models[['TVAR1']] <- data_train %>% 
  select(-date) %>% 
  ts() %>% 
  tsDyn::TVAR(lag = 1, include = 'const', model = 'TAR',
            mTh = which(names(select(data, -date)) == 'RP'))

fitted(models[['TVAR1']])

#TVAR2
models[['TVAR2']] <- data_train %>% 
  select(-date) %>% 
  ts() %>% 
  tsDyn::TVAR(lag = 1, include = 'const', model = 'TAR',
              mTh = which(names(select(data, -date)) == 'RP'),
              nthresh = 2)

fitted(models[['TVAR2']])



# Out of sample forecast

forecast <- bind_rows(
forecast(models[['AR1']], 12) %>% 
  as_tibble() %>% 
  select(value = `Point Forecast`) %>% 
  rownames_to_column() %>% 
  mutate(model = 'AR1'),

forecast(models[['ARiMA']], 12) %>% 
  as_tibble() %>% 
  select(value = `Point Forecast`) %>% 
  rownames_to_column() %>% 
  mutate(model = 'ARiMA'),

forecast(models[['VAR']], 12)$forecast$HP %>% 
  as_tibble() %>% 
  select(value = `Point Forecast`) %>% 
  rownames_to_column() %>% 
  mutate(model = 'VAR'),

predict(models[['TVAR1']], n.ahead = 12) %>% 
  as_tibble() %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'TVAR1'),

predict(models[['TVAR2']], n.ahead = 12) %>% 
  as_tibble() %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'TVAR2'),

data_test %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'Actual')

) %>% 
  mutate(rowname = as.numeric(rowname))

forecast  %>% 
  ggplot(aes(x = rowname, y = value, color= model, group = model)) +
  geom_line(linewidth = 1)


#Fitted
fitted <- bind_rows(
fitted(models[['AR1']]) %>% 
  as_tibble() %>% 
  rename(value = x) %>% 
  rownames_to_column() %>% 
  mutate(model = 'AR1'),

fitted(models[['ARiMA']]) %>% 
  as_tibble() %>% 
  rename(value = x) %>% 
  rownames_to_column() %>% 
  mutate(model = 'ARiMA'),

bind_rows(
rep(NA, models[['VAR']]$p) %>% 
  as_tibble() %>% 
  mutate(value = as.numeric(value)),
fitted(models[['VAR']]) %>% 
  as_tibble() %>% 
  select(value = HP)
) %>% 
  rownames_to_column() %>% 
  mutate(model = 'VAR'),

bind_rows(
  rep(NA, models[['TVAR1']]$lag) %>% 
    as_tibble() %>% 
    mutate(value = as.numeric(value)),
  fitted(models[['TVAR1']]) %>% 
    as_tibble() %>% 
    select(value = HP)
) %>% 
  rownames_to_column() %>% 
  mutate(model = 'TVAR1'),

bind_rows(
  rep(NA, models[['TVAR2']]$lag) %>% 
    as_tibble() %>% 
    mutate(value = as.numeric(value)),
  fitted(models[['TVAR2']]) %>% 
    as_tibble() %>% 
    select(value = HP)
) %>% 
  rownames_to_column() %>% 
  mutate(model = 'TVAR2'),

data_train %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'Actual')

) %>% 
  mutate(rowname = as.numeric(rowname)) 

fitted %>% 
  ggplot(aes(x = rowname, y = value, color= model, group = model)) +
  geom_line(linewidth = 1)



