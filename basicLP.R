data <- ir %>% 
  #filter(indicator == 'FIMM_PA') %>% 
  #unnest(data) %>% 
  #ungroup() %>% 
  select(-indicator, -value, -contains('AR')) %>% 
  inner_join(hprice %>% 
               select(-value), 
             by = c('ccode2', 'date'))

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

data %>% 
  group_by(ind) %>% 
  mutate(ind = str_remove_all(ind, 'q') %>% as.numeric()) %>% 
  nest() %>% 
  mutate(lpirf = map(data, ~gen_lpirf(.x))) %>% 
  select(lpirf) %>% 
  unnest(lpirf) %>% 
  ggplot(aes(x = t, y = mean, color = ind)) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = ind), alpha = .2, linetype = 'dashed') +
  geom_hline(aes(yintercept = 0), color = 'red') +
  facet_wrap(~ind, scales = 'free') +
  theme_minimal()

