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
models[['BVAR']] <- BVAR::bvar(data = data_train,
                               lags = 12,
                               n_draw = 60000,
                               n_burn = 10000,
                               priors = bv_priors())



#TVAR
models[['TVPVAR']] <- shrinkTVPVAR::shrinkTVPVAR(y = data_train %>% 
                                                   select(-date, -RP, -R, 
                                                          -EMP, -POP, -C) %>% 
                                                   ts(),
                                                 p = 12,
                                                 mod_type = 'triple',
                                                 niter = 60000,
                                                 nburn = 10000)

# Out of sample forecast

fct_tvp_hp <- function(tvpvar){
  fct <- shrinkTVPVAR::forecast_shrinkTVPVAR(tvpvar, 12)$HP_forc$y_pred
  out <- matrix(nrow = 12, ncol = 1)
  for(i in 1:ncol(fct)){
    out[i,] <- fct[,i] %>% 
      median()
  }
  colnames(out) = 'value'
  out <- out %>% 
    as_tibble() %>% 
    rownames_to_column() %>% 
    mutate(model = 'TVPVAR')
  out
}


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

bvars::forecast(models[['BVAR']], 12)$forecast %>% 
  as_tibble() %>% 
  drop_na() %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'BVAR'),

fct_tvp_hp(models[['TVPVAR']]),

data_test %>% 
  select(value = HP) %>% 
  rownames_to_column() %>% 
  mutate(model = 'Actual')

) %>% 
  mutate(rowname = as.numeric(rowname))

forecast  %>% 
  ggplot(aes(x = rowname, y = value, color= model, group = model)) +
  geom_line(linewidth = 1)



