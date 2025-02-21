
data_prep <- data_shocks %>% 
  group_by(ccode2) %>% 
  mutate(across(c(-date, -contains("shock")), ~.x-mean(.x, na.rm = T))) %>% 
  #mutate(lvalue = RPPI) %>% 
  #mutate(lvalue = RPPI - lag(RPPI, 4)) %>% 
  mutate(lvalue = RPPI - lag(RPPI)) %>%
  drop_na() %>% 
  mutate(ind = case_when(#lvalue <= quantile(lvalue, .1) ~ 'q1',
                         lvalue <= quantile(lvalue, .2) ~ 'q2',
                         #lvalue <= quantile(lvalue, .3) ~ 'q3',
                         lvalue <= quantile(lvalue, .4) ~ 'q4',
                         #lvalue <= quantile(lvalue, .5) ~ 'q5',
                         lvalue <= quantile(lvalue, .6) ~ 'q6',
                         #lvalue <= quantile(lvalue, .7) ~ 'q7',
                         lvalue <= quantile(lvalue, .8) ~ 'q8',
                         #lvalue <= quantile(lvalue, .9) ~ 'q9',
                         lvalue <= quantile(lvalue, 1) ~ 'q10')) %>% 
  select(-lvalue) %>% 
  gather(key = var, value = value,
         -ccode2, -date, -ind, -contains('shock')) %>% 
  group_by(ccode2, var) %>% 
  nest() %>% 
  mutate(data = map(data, ~.x %>% 
                      #Create RHS
                      mutate(map_dfc(seq(4), ~ lag(value, n = .x)) %>%
                               set_names(paste('lag', seq(4),sep = ''))) %>% 
                      #Create RHS
                      mutate(map_dfc(seq(4), ~ lag(shock, n = .x)) %>%
                               set_names(paste('shock_lag', seq(4),sep = ''))) %>%
                      #Create LHS
                      mutate(map_dfc(seq(21), ~ lead(value, n = .x)) %>%
                               set_names(paste('lead', seq(21),sep = ''))) %>% 
                      select(-value)
                      )) 


gen_lpirf <- function(data){
  
  require(tidyverse)
  require(stringr)
  require(lubridate)
  require(rebus)
  require(sandwich)
  
  
  x <- data %>% 
    #drop_na() %>% 
    ungroup() %>% 
    select(contains('shock') | contains('lag'))
  x_names <- names(x)
  shocks_names <- x_names[!str_detect(x_names, 'lag')]
  x <- as.matrix(x)
  
  y <- data %>% 
    #drop_na() %>% 
    ungroup() %>% 
    select(contains('lead'))
  y_names <- names(y)
  y <- as.matrix(y)
  
  
  m <- round(0.75*sqrt(nrow(x)))
  h <-  names(data)[str_detect(names(data), 'lead')] %>% 
    str_extract_all(one_or_more(DGT)) %>% 
    as.matrix() %>% 
    as.numeric() %>% 
    max()
  
  
  
  #Estimate the IRF for one variable
  out <- NULL
  pb = txtProgressBar(min = 0, max = h, initial = 0) 
  for(i in 1:h){
    
    setTxtProgressBar(pb,i)
    
    model <- lm(y[,i] ~ -1 + x)
    irf <- model$coefficients[!str_detect(names(model$coefficients), or('lag', 'Intercept', 'dist'))]
    names(irf) <- shocks_names
    sum <- summary(model)
    se <- sum$coefficients
    se <- se[rownames(se)[!str_detect(rownames(se),or('lag', 'Intercept', 'dist'))], 
             colnames(se)[colnames(se) == 'Std. Error']]
    
    #nw <- NeweyWest(model, 
    #                lag = m - 1, 
    #                prewhite = F, 
    #                adjust = T)
    #se <- sqrt(diag(nw))[1:(ncol(x))]
    #se <- se[!str_detect(names(se), or('lag', 'Intercept'))]
    names(se) <- shocks_names
    irf_ub <- irf + se
    names(irf_ub) <- paste(names(irf_ub), '_ub', sep = '')
    irf_lb <- irf - se
    names(irf_lb) <- paste(names(irf_lb), '_lb', sep = '')
    out[[i]] <- as_tibble(t(c(irf, irf_ub, irf_lb))) %>% 
      gather(key = key,
             value = value) %>% 
      mutate(key2 = case_when(str_detect(key, 'ub') ~ 'ub',
                              str_detect(key, 'lb') ~ 'lb',
                              TRUE ~ 'mean'),
             key = str_remove_all(key, key2) %>% 
               str_remove_all('_' %R% END)) %>% 
      spread(key2, value) %>% 
      rename(shock = key) %>% 
      mutate(t = i)
  }
  close(pb)
  out <- bind_rows(out) %>% 
    mutate(t = t - 1) 
  
  gc()
  
  out
}


lpirfs <- data_prep %>% 
  unnest(data) %>% 
  ungroup() %>% 
  group_by(var, ind) %>% 
  select(-date) %>% 
  nest() %>% 
  mutate(data1 = map(data, ~.x %>% 
                       select(-contains('shock2')))) %>%
  mutate(lpirf1 = map(data1, ~gen_lpirf(.x))) 


lpirfs %>%
  mutate(ind = str_remove_all(ind, 'q') %>% as.numeric()) %>% 
  select(var, ind, lpirf1) %>% 
  unnest(lpirf1) %>%
  filter(var == 'RPPI') %>% 
  ggplot(aes(x = t, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = ind), alpha = .2, linetype = 'dashed') +
  geom_hline(aes(yintercept = 0), color = 'red') +
  facet_wrap(~ind, scales = 'free') +
  theme_minimal() +
  theme(legend.position = 'none')










