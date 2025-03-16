x1 <- data %>% 
  select(HP) %>% 
  ts(frequency = 4,
     start = c(data$date %>% min %>% year,
               data$date %>% min %>% quarter)) %>% 
  BCDating::BBQ(mincycle = 1,
                minphase = 1) 


BCDating::summary(x1)

BCDating::avgts(data %>% 
                  select(HP) %>% 
                  ts(frequency = 4,
                     start = c(data$date %>% min %>% year,
                               data$date %>% min %>% quarter)),
                x1) %>% 
  print() %>% 
  plot()

BCDating::plot(x1, data %>% 
                 select(HP) %>% 
                 ts(frequency = 4,
                    start = c(data$date %>% min %>% year,
                              data$date %>% min %>% quarter)))

