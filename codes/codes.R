library(fredr)
library(tidyverse)
library(vars)
library(tsDyn)
library(urca)
library(lubridate)
library(remotes)
install_github("angusmoore/tvarGIRF")
library(tvarGIRF)


#import data
#set api key
fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")

#download data
params <- list(
  series_id = c("QUSR628BIS", "MORTGAGE30US", "DFF"),
  frequency = "q",
  observation_start = as.Date("1950-01-01")
)


import  <- pmap_dfr(
  .l = params,
  .f = ~ fredr(series_id = .x, frequency = .y)
) %>%
  dplyr::select(date, series_id, value) %>%
  spread(key = series_id, value = value) %>%
  drop_na() %>% rename(ffr = DFF,
                       m30 = MORTGAGE30US,
                       hprice = QUSR628BIS) %>%
  drop_na() %>% 
  filter(date <= as.Date("2021-12-01")) %>% 
  left_join(
fredr("GDPC1",
             frequency = "q") %>% 
  dplyr::select(date, gdp = value) %>% 
  mutate(year = year(date)) %>% 
  group_by(year) %>% 
  mutate(mean = mean(gdp),
         mean_2010 = ifelse(year == 2010, mean, NA)) %>% 
  ungroup() %>% 
  mutate(mean_2010 = mean(mean_2010, na.rm = TRUE),
         gdp = 100*gdp/mean_2010) %>% 
  dplyr::select(date, gdp), by = "date"
) %>% 
left_join(
  fredr("RELACBW027SBOG",
                frequency = "q") %>% 
            dplyr::select(date, hloan = value) %>% 
            mutate(year = year(date)) %>% 
            group_by(year) %>% 
            mutate(mean = mean(hloan),
                   mean_2010 = ifelse(year == 2010, mean, NA)) %>% 
            ungroup() %>% 
            mutate(mean_2010 = mean(mean_2010, na.rm = TRUE),
                   hloan = 100*hloan/mean_2010) %>% 
            dplyr::select(date, hloan),
  by = "date") %>% 
  drop_na()



import %>% gather(key = "variable", value = "value", ffr, m30, hloan, gdp, hprice) %>%
  mutate(variable = case_when(variable == "ffr" ~ "Federal Funds Rate",
                              variable == "gdp" ~ "GDP",
                              variable == "hloan" ~ "Real Estate Loans",
                              variable == "hprice" ~ "Real Residential Property Prices",
                              variable == "m30" ~ "30-Year Fixed Rate Mortgage Average")) %>%
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~variable, scales = "free") + theme_minimal() +
  labs(x = "",
       y = "",
       caption = "The real GDP and the volume of real estate loans are normalized by their 2010 average") 



data <- import %>% 
  dplyr::select(hprice, ffr, m30, hloan, gdp) %>% ts()


VARselect(data)$selection %>% as_tibble() %>% 
  mutate(type = names(VARselect(data)$selection)) %>% 
  mutate(type = str_replace_all(type, "\\(n\\)", "")) %>% 
  dplyr::select(type, value)


johansen <- ca.jo(data, type = "trace", K = 2, spec = "longrun", ecdet = "none") %>% 
  summary()


coint_tab <- johansen@cval %>% 
  as.tibble() %>% 
  mutate(test = johansen@teststat) %>% 
  dplyr::select(test, everything()) %>% 
  as.matrix()
rownames(coint_tab) <- rownames(johansen@cval)
for(i in 1:nrow(coint_tab)){
  for(j in 1:ncol(coint_tab)){
    coint_tab[i,j] <- round(coint_tab[i,j], digits = 2)
  }
    }
coint_tab <- coint_tab %>% stargazer::stargazer(title = "Trace statistics from Johansen cointegration test",
                                   style = "aer",
                                   notes = "With two lags included as suggested by the Schwarz-Bayesian Information Criterion, and constant")
coint_tab[4:length(coint_tab)]


vec <- VECM(data, lag  =2, r = 1, estim = "2OLS", include = "none") %>% 
  summary()

vec$model.specific$coint[2:length(vec$model.specific$coint)]


TVECM.SeoTest(data, lag = 2, beta = -1*c(-1.664798, -0.7307156, 0.1000535, -1.168652),
              nboot = 10, plot = TRUE)





data <- import %>% 
  dplyr::select(ffr, m30, hloan, gdp, hprice) %>% 
  mutate(hloan = hloan - lag(hloan),
         gdp = gdp - lag(gdp),
         hprice = hprice - lag(hprice)) %>% drop_na() %>% ts()


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
#e <- resid(var)
#cov_mat <- t(e) %*% e
#chol <- chol(cov_mat)
#ffrshock <- chol[1,] / chol[1,1]

source("SVAR2.R")
struct_shock <- function(var, data, amat, bmat, estmethod){
  x <- SVAR2(var = var, data = data,
             Amat = amat,
             Bmat = bmat, 
             estmethod = estmethod)
  
  shockmat <- x$A %*% x$B
  
  shockmat <- shockmat*-1
  diag(shockmat) <- diag(shockmat)*-1
  
  shockmat <- shockmat / diag(shockmat)
  
  message(paste("Minimum number of restrictions in the A matrix must be ", ncol(data)*(ncol(data)-1)/2, sep = ""))
  shockmat
}

amat <- diag(ncol(data))
amat[lower.tri(amat)] <- NA
diag(amat) <- NA
#FFR shocks
amat[4:5, 1] <- 0
#M30 shocks
amat[4, 2] <- 0
#HLOAN shocks
amat[4, 3] <- 0
#GDP shocks
amat[5,4] <- 0
amat[1,4] <- NA

bmat <- diag(ncol(data))
diag(bmat) <- NA

ffrshock <- struct_shock(var = var,
                         data = data,
                         amat = amat,
                         bmat = bmat,
                         estmethod = "scoring")

ffrshock <- ffrshock[,1]


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
       nahead = 40,
       coefmat = varcoef,
       main = "Impulse responses of FFR shock, simple SVAR(1)")


#tesztelj?k a tvar helyess?g?t
par_orig <- c(5.1, 4.1, 4.1, 2.1)
par_large <- c(2,2,2,2)
par(mar = par_large)
TVAR.LRtest(data, lag = 1, mTh = 5, plot = TRUE, nboot = 1000)
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
#e <- resid(tvar)
#cov_mat <- t(e) %*% e
#chol <- chol(cov_mat)
#ffrshock <- chol[1,] / chol[1,1]

ffrshock <- struct_shock(var = tvar,
                         data = data,
                         amat = amat,
                         bmat = bmat,
                         estmethod = "scoring")

ffrshock <- ffrshock[,1]


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
       nahead = 40,
       coefmat = highcoef,
       main = "Impulse responses from FFR shock, TVAR(1) - high regime")

irfgen(shock = ffrshock,
       nahead = 40,
       coefmat = lowcoef,
       main = "Impulse responses from FFR shock, TVAR(1) - low regime")

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
  geom_line(linewidth = 1) +
  facet_wrap(~variable, scales = "free") + theme_minimal() +
  scale_color_manual(breaks = c("High regime", "Low regime"),
                     values = c("#E9002D", "#00B000")) +
  labs(x = "",
       y = "") +
  theme(legend.title = element_blank()) +
  geom_rect(data = fredr(series_id = "USREC",
                         frequency = "m") %>% 
              dplyr::select(date, recession = value) %>% 
              mutate(diff = recession - lag(recession)) %>% 
              filter(!is.na(diff)) %>% 
              mutate(recession_start = ifelse(diff == 1, as.character(date), NA),
                     recession_end = ifelse(diff == -1, as.character(date), NA)) %>% 
              filter(!is.na(recession_start) | !is.na(recession_end)) %>% 
              mutate(recession_end = ifelse(!is.na(recession_start), lead(recession_end), NA)) %>%
              filter(!is.na(recession_start)) %>% 
              mutate(across(.cols = c(recession_start, recession_end), .fns = ~as.Date(.x))) %>% 
              dplyr::select(recession_start, recession_end) %>% 
              filter(recession_start >= min(import$date),
                     recession_start <= max(import$date)),
            inherit.aes = F,
            aes(xmin = recession_start, xmax = recession_end, ymin = -Inf, ymax = Inf), 
            fill = "grey50", alpha = 0.5)



GIRF <- function(tvar, horizon, reps, hists, shock, seed,
                 aggregate = c("mean", "median"), prob = NA, group = FALSE){
set.seed(seed)
irfs <- NULL

regime_track <- bind_cols(1:tvar$T,
          tvar$model.specific$regime) %>% 
  rename(obs = ...1,
         regime = ...2) %>% 
  filter(!is.na(regime))

resid <- tvar$residuals

resids <- regime_track %>% bind_cols(resid %>% as_tibble())


start <- (!is.na(tvar$model.specific$regime)) %>% which() %>% min()

for(k in start:nrow(data)){
  gc()
  obs_num <- k
  regime <- tvar$model.specific$regime[k]
  tvec <- 1:horizon
  
  shockmat <- matrix(0, nrow = horizon, ncol = tvar$k)
  shockmat[1, ] <- shock
  
  resmat <- matrix(0, nrow = horizon, ncol = tvar$k)
  
  r <- regime_track %>% filter(obs == k) %>% pull(regime)
  res <- resids %>% filter(regime == r) %>% dplyr::select(-obs, - regime)
    
Y2 <- matrix(data = 0, nrow = horizon, ncol = tvar$k)
for(j in 1:hists){
Y <- matrix(data = 0, nrow = horizon, ncol = tvar$k)
for(i in 1:reps){
  res <- res[sample(nrow(res), 1), ] %>% as.matrix()
  resmat[1, ] <- res
  Y_null <- TVAR.sim(B = tvar$coeffmat,
                     Thresh = tvar$model.specific$Thresh,
                     nthresh = tvar$model.specific$nthresh,
                     n = horizon,
                     lag = tvar$lag,
                     thDelay = tvar$model.specific$thDelay,
                     mTh = 5,
                     starting = data[(k:k+tvar$lag-1),] %>% t(),
                     innov = resmat)

Y_shock <- TVAR.sim(B = tvar$coeffmat,
                    Thresh = tvar$model.specific$Thresh,
                    nthresh = tvar$model.specific$nthresh,
                    n = horizon,
                    lag = tvar$lag,
                    thDelay = tvar$model.specific$thDelay,
                    mTh = 5,
                    starting = data[(k:k+tvar$lag-1),] %>% t(),
                    innov = shockmat)


Y <- Y + (Y_shock - Y_null)/reps
}
Y2 <- Y2 + Y/hists
}
  irfdata <- cbind(tvec, Y2)
  colnames(irfdata) <- c("t", colnames(data))
  irfdata <- as_tibble(irfdata) %>% list()

  
  irfs[[k]] <- tibble(obs_num, regime, irfdata)
}

irfs <- bind_rows(irfs) %>% 
  filter(!is.na(regime)) %>% 
  unnest(irfdata)

variables <- colnames(irfs)[4:length(colnames(irfs))]

if(group == TRUE){
  irfs <- irfs %>% group_by(regime, t)
} else if(group == FALSE) {
  irfs <- irfs %>% group_by(t)
}

if(aggregate == "mean"){
  irf_c <- irfs %>% 
    summarize(across(.cols = all_of(variables), .fns = ~mean(.x))) %>% 
    gather(key = "variable", value = "response",
           all_of(variables))
} else if(aggregate == "median"){
  irf_c <- irfs %>% 
    summarize(across(.cols = all_of(variables), .fns = ~median(.x))) %>% 
    gather(key = "variable", value = "response",
           all_of(variables))
}




if(!is.na(prob)){
  irf_u <- irfs %>% 
    summarize(across(.cols = all_of(variables), .fns = ~quantile(.x, 1-prob))) %>% 
    gather(key = "variable", value = "upper",
           variables)
  irf_l <- irfs %>% 
    summarize(across(.cols = all_of(variables), .fns = ~quantile(.x, prob))) %>% 
    gather(key = "variable", value = "lower",
           variables)
  
  irfs <- irf_c %>% 
    left_join(irf_l) %>% 
    left_join(irf_u)
  
} else if(is.na(prob)){
  irfs <- irf_c
}

irfs

}

irfs1 <- GIRF(tvar = tvar,
     horizon = 40,
     reps = 100,
     hists = 100,
     shock = ffrshock,
     seed = 1234,
     aggregate = "median",
     prob = 0.16)


p1 <- irfs1 %>% 
  mutate(variable = case_when(variable == "ffr" ~ "FFR",
                              variable == "m30" ~ "Mortgage rate",
                              variable == "hloan" ~ "Housing loan",
                              variable == "gdp" ~ "GDP",
                              variable == "hprice" ~ "House Price Index"),
         variable = factor(variable, levels = c("FFR", "Mortgage rate", "Housing loan", "GDP", "House Price Index"))) %>%
  ggplot(aes(x = t, y = response)) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red")+
#  geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5) +
  facet_wrap(~variable, scales = "free")

irfs2 <- GIRF(tvar = tvar,
              horizon = 40,
              reps = 100,
              hists = 100,
              shock = ffrshock,
              seed = 1234,
              aggregate = "median",
              prob = 0.16,
              group = TRUE)

p2 <- irfs2 %>% 
  mutate(regime = case_when(regime == 1 ~ "Low regime",
                               regime == 2 ~ "High regime")) %>% 
  mutate(variable = case_when(variable == "ffr" ~ "FFR",
                              variable == "m30" ~ "Mortgage rate",
                              variable == "hloan" ~ "Housing loan",
                              variable == "gdp" ~ "GDP",
                              variable == "hprice" ~ "House Price Index"),
         variable = factor(variable, levels = c("FFR", "Mortgage rate", "Housing loan", "GDP", "House Price Index"))) %>% 
  ggplot(aes(x = t, y = response, color = as.factor(regime))) +
  geom_line() +
  geom_hline(yintercept = 0, color = "red")+
  facet_wrap(~variable, scales = "free")+
  scale_color_manual(breaks = c("High regime", "Low regime"),
                     values = c("#E9002D", "#00B000")) +
  labs(x = "",
       y = "")+
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11),
        legend.title = element_blank())



p1 +
  geom_smooth(method = "gam", se = F)

p2 +
  geom_smooth(method = "gam", se = F)

irfs1 %>% 
  saveRDS("GIRF_1.rds")
irfs2 %>% 
  saveRDS("GIRF_2.rds")

p1 %>% ggsave("girfs1.png")
p2 %>% ggsave("girfs2.png")



