#packages
suppressPackageStartupMessages({
library(readxl)
library(tidyverse)
library(tseries)
library(forecast)
library(TSA)
library(aTSA)
library(fpp2)
library(lmtest)
library(vars)
library(mFilter)
library(ggplot2)
library(tsDyn)
library(fredr)
library(lubridate)
library(cowplot)
library(gridExtra)
library(patchwork)
library(data.table)
library(tsibble)
library(xts)
library(zoo)
library(broom)
library(urca)
})



folder <- "C:/Users/horva/Documents/Egyetem/PhD/Study material/Semester 2/Time series/Empirical project"

#set api key
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

#download data
params <- list(
  series_id = c("QUSR628BIS", "RELACBW027SBOG", "MORTGAGE30US", "DFF", "GDPC1"),
  frequency = "q",
  observation_start = as.Date("1950-01-01")
)


import  <- pmap_dfr(
  .l = params,
  .f = ~ fredr(series_id = .x, frequency = .y)
) %>%
  dplyr::select(date, series_id, value) %>%
  spread(key = series_id, value = value) %>%
  drop_na() %>% as_tsibble() %>% rename(ffr = DFF,
                                        m30 = MORTGAGE30US,
                                        hloan = RELACBW027SBOG,
                                        gdp = GDPC1,
                                        hprice = QUSR628BIS) %>%
  drop_na()

import %>% as_tibble() %>% gather(key = "variable", value = "value", ffr, m30, hloan, gdp, hprice) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "gdp" ~ "GDP",
                              variable == "hloan" ~ "Real Estate Loans",
                              variable == "hprice" ~ "Real Residential Property Prices",
                              variable == "m30" ~ "30-Year Fixed Rate Mortgage Average")) %>%
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free") + theme_minimal() +
  labs(x = "",
       y = "") 
  ggsave(file.path(folder, "dataplot.png"),
         device = "png")


data <- import %>% 
  mutate(hloan = log(hloan),
         gdp = log(gdp),
         hprice = log(hprice),
         hloan = hloan - lag(hloan),
         gdp = gdp - lag(gdp),
         hprice = hprice - lag(hprice)) %>%
  drop_na() %>%
  relocate(ffr, .after = date) %>%
  relocate(m30, .after = ffr) %>%
  relocate(hloan, .after = m30) %>%
  relocate(gdp, .after = hloan) %>%
  relocate(hprice, .after = gdp)


#plotting log-differenced gdp, house price index and housing loand volume
ggplot(data) +
  geom_line(aes(x = date, y = gdp), color = "blue") +
  geom_line(aes(x = date, y = hprice), color = "darkgreen") +
  geom_line(aes(x = date, y = hloan), color = "red") +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11)) +
  labs(y = "",
       x = "Date",
       title = "GDP (blue), House Price Index (green) and Housing loan volume (red) for the US, log-difference")

#cross correlations
cor(data$gdp, data$hprice)
cor(data$gdp, data$hloan)
cor(data$hloan, data$hprice)

data <- data[, 2:6] %>% ts()




#fitting simple var model
VARselect(data)

var <- VAR(data, p = 1, type = "const")

varcoef <- bind_rows(
var$varresult$ffr$coefficients,
var$varresult$m30$coefficients,
var$varresult$hloan$coefficients,
var$varresult$gdp$coefficients,
var$varresult$hprice$coefficients) %>%
  as.matrix()

#identify structural shocks
e <- resid(var)
cov_mat <- t(e) %*% e
chol <- chol(cov_mat)

print(chol) %>% t()

ffrshock <- chol[1,] / chol[1,1]

irfgen <- function(shock, nahead, coefmat, main){
  
  irf<- matrix(nrow = length(shock), ncol = nahead+1)
  irf[,1] <- shock %>% as.matrix()
  t <- matrix(ncol = 1, nrow = nahead+1)
  
  for(j in 1:ncol(irf)){
    t[j, 1] <- j-1
  }
  
  for(j in 2:ncol(irf)){
    for(i in 1:nrow(irf)){
      irf[i,j] <- irf[1,j-1]*coefmat[i,1]+
        irf[2,j-1]*coefmat[i,2]+
        irf[3,j-1]*coefmat[i,3]+
        irf[4,j-1]*coefmat[i,4]+
        irf[5,j-1]*coefmat[i,5]
    }
  }
  
  irf <- t(irf)
  
  colnames(irf) <- names(shock)
  irf <- bind_cols(t, irf) %>%
    as_tibble() %>%
    rename(t = ...1)
  
  irf <- gather(irf, key = "variable", value = "response", ffr, m30, hloan, gdp, hprice) %>%
    mutate(variable = case_when(variable == "ffr" ~ "FFR",
                                variable == "m30" ~ "Mortgage rate",
                                variable == "hloan" ~ "Housing loan",
                                variable == "gdp" ~ "GDP",
                                variable == "hprice" ~ "House Price Index"),
           variable = factor(variable, levels = c("FFR", "Mortgage rate", "Housing loan", "GDP", "House Price Index")))
  
  ggplot(irf, aes(x = t, y = response)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "red")+
    facet_wrap(~variable, scales = "free") +
    labs(x = "",
         y = "",
         title = main)+
    theme(plot.title = element_text(size = 11, hjust=0.5),
          axis.title.y = element_text(size=11))
  
}

irfgen(shock = ffrshock,
       nahead = 20,
       coefmat = varcoef,
       main = "Impulse responses of FFR shock, simple SVAR(1)")
ggsave(file.path(folder, "svarirf.png"),
       device = "png")

#teszteljük a tvar helyességét
par_orig <- c(5.1, 4.1, 4.1, 2.1)
par_large <- c(2,2,2,2)
par(mar = par_large)
png(file.path(folder, "tvartest.png"))
TVAR.LRtest(data, lag = 1, mTh = 5, plot = TRUE, nboot = 10000)
dev.off()
par(mar = par_orig)

#tvar with optimal threshold value - possibly not very useful
tvar <- TVAR(data, lag = 1, nthresh = 1, mTh = 5, model = "TAR", max.iter = 1000, trim = 0)  
print(tvar)
summary(tvar)
plot(tvar)

#tvar with threshold value = 0 - a compromised solution
tvar <- TVAR(data, lag = 1, nthresh = 1, mTh = 5, model = "TAR", max.iter = 1000, gamma = 0)  
print(tvar)
summary(tvar)
plot(tvar)

#identify structural shocks
e <- resid(tvar)
cov_mat <- t(e) %*% e
chol <- chol(cov_mat)

print(chol) %>% t()

ffrshock <- chol[1,] / chol[1,1]

highcoef <- tvar$coefficients$Bup

lowcoef <- tvar$coefficients$Bdown

irfgen <- function(shock, nahead, coefmat, main){
  
  irf<- matrix(nrow = length(shock), ncol = nahead+1)
  irf[,1] <- shock %>% as.matrix()
  t <- matrix(ncol = 1, nrow = nahead+1)
  
  for(j in 1:ncol(irf)){
    t[j, 1] <- j-1
  }
  
  for(j in 2:ncol(irf)){
    for(i in 1:nrow(irf)){
      irf[i,j] <- irf[1,j-1]*coefmat[i,2]+
        irf[2,j-1]*coefmat[i,3]+
        irf[3,j-1]*coefmat[i,4]+
        irf[4,j-1]*coefmat[i,5]+
        irf[5,j-1]*coefmat[i,6]
    }
  }
  
  irf <- t(irf)
  
  colnames(irf) <- names(shock)
  irf <- bind_cols(t, irf) %>%
    as_tibble() %>%
    rename(t = ...1)
  
  irf <- gather(irf, key = "variable", value = "response", ffr, m30, hloan, gdp, hprice) %>%
    mutate(variable = case_when(variable == "ffr" ~ "FFR",
                                variable == "m30" ~ "Mortgage rate",
                                variable == "hloan" ~ "Housing loan",
                                variable == "gdp" ~ "GDP",
                                variable == "hprice" ~ "House Price Index"),
           variable = factor(variable, levels = c("FFR", "Mortgage rate", "Housing loan", "GDP", "House Price Index")))
  
  ggplot(irf, aes(x = t, y = response)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "red")+
    facet_wrap(~variable, scales = "free")+
    labs(x = "",
         y = "",
         title = main)+
    theme(plot.title = element_text(size = 11, hjust=0.5),
          axis.title.y = element_text(size=11))
  
  
}


irfgen(shock = ffrshock,
                   nahead = 20,
                   coefmat = highcoef,
       main = "Impulse responses from FFR shock, TVAR(1) - high regime")
ggsave(file.path(folder, "tvarhigh.png"),
       device = "png")

irfgen(shock = ffrshock,
       nahead = 20,
       coefmat = lowcoef,
       main = "Impulse responses from FFR shock, TVAR(1) - low regime")
ggsave(file.path(folder, "tvarlow.png"),
       device = "png")

indicator <- tvar$model.specific$regime

import %>% filter(date >= as.Date("1973-04-01")) %>%
  bind_cols(indicator) %>%
  drop_na() %>%
  rename(indicator = ...7) %>%
  as_tibble() %>% gather(key = "variable", value = "value", ffr, m30, hloan, gdp, hprice) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "gdp" ~ "GDP",
                              variable == "hloan" ~ "Real Estate Loans",
                              variable == "hprice" ~ "Real Residential Property Prices",
                              variable == "m30" ~ "30-Year Fixed Rate Mortgage Average")) %>%
  mutate(indicator = case_when(indicator == 1 ~ "Low regime",
                               indicator == 2 ~ "High regime")) %>%
  ggplot(aes(x = date, y = value, color = indicator, group = 1)) +
  geom_line(size = 1) +
  facet_wrap(~variable, scales = "free") + theme_minimal() +
  labs(x = "",
       y = "") +
  theme(legend.title = element_blank()) 
ggsave(file.path(folder, "indicator.png"),
       device = "png")

