library(fredr)
library(tidyverse)
library(lubridate)
library(seasonal)
library(rlang)
library(forecast)
library(gridExtra)
library(tsDyn)
library(parallel)
library(meboot)

source('codes/data.R')
source('codes/tvarestfuns.R')
data <- data()

ts_orig <- data %>% 
  select(Y, R, TS, CS, MR, HL, CGG, P, HP, RP) %>% 
  ts(freq = 4, 
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min))

createfactors <- function(data){
allvars <- data %>% 
  select(Y, R, TS, CS, MR, HL, CGG, P, HP, RP) %>% 
  ts(freq = 4, 
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min)) %>% 
  colnames()

RPCGG <- data %>% 
  select(RP, CGG) %>% 
  mutate(across(c(RP, CGG), ~scale(.x)[,1])) %>% 
  FactoMineR::PCA(graph = F) 

loadings1 <- RPCGG$var$coord

RPCGG <- RPCGG$ind$coord[,1] %>% 
  ts(freq = 4, 
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min))


X <- data %>% 
  select(-CGG,-R, -RP, -date) %>% 
  mutate(across(everything(), ~scale(.x)[,1])) %>% 
  FactoMineR::PCA(graph = F) 


loadings2 <- X$var$coord[,1:3]

X <- X$ind$coord[,1:3] %>% 
  ts(freq = 4, 
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min))

R <- data %>% 
  select(R) %>% 
  ts(freq = 4, 
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min))

ts <- ts.union(R = R, RPCGG, X)
colnames(ts) <- c('R', 'RPCGG', 'F1', 'F2', 'F3')


loadings_full <- matrix(0,
                        nrow = length(allvars),
                        ncol = ncol(ts))
colnames(loadings_full) <- colnames(ts)
rownames(loadings_full) <- allvars
loadings_full['R', 'R'] <- 1
loadings_full[rownames(loadings1), 'RPCGG'] <- loadings1[,1]
loadings_full[rownames(loadings2), colnames(loadings_full) != c('R', 'RPCGG')] <- loadings2
loadings_full

out <- list(ts = ts,
            loadings_full = loadings_full)

out
}
  

data %>% 
  select(-CGG,-R, -RP, -date) %>% 
  mutate(across(everything(), ~scale(.x)[,1])) %>% 
  FactoMineR::PCA(graph = F) %>% 
  summary()

(data %>% 
  select(-CGG,-R, -RP, -date) %>% 
  mutate(across(everything(), ~scale(.x)[,1])) %>% 
  FactoMineR::PCA(graph = F))$var$coord 

(data %>% 
  select(-CGG,-R, -RP, -date) %>% 
  mutate(across(everything(), ~scale(.x)[,1])) %>% 
  FactoMineR::PCA(graph = F))$var$contrib 

weights <- PCA_weights(data)

makePCAplot(data, weights)

#LRtests(ts, mTh = 2)

ts <- createfactors(data)$ts
loadings_full <- createfactors(data)$loadings_full

thvar_compar(ts, data, weights)

tvar <- tsDyn::TVAR(ts, 
                    lag = 4,
                    nthresh = 1,
                    mTh = 2)


gen_sr <- function(loadings_full){
sign_irf <- loadings_full/abs(10^-100 + loadings_full)
for(i in 1:nrow(sign_irf)){
  for(j in 1:ncol(sign_irf)){
    sign_irf[i,j] <- ifelse(sign_irf[i,j] == 0, NA, sign_irf[i,j])
  }
}
sign_irf[c('R', 'MR', 'CS', 'RP'),'R'] <- 1
sign_irf[c('P', 'HP', 'HL', 'TS', 'CGG'),'R'] <- -1
#message('Applying the following set of sign restrictions:')
#print(sign_irf)
colnames(sign_irf) <- NULL
rownames(sign_irf) <- NULL
#message('To modify the resttrictions, please edit the function.')
#message('Storing sign restrictions without dimnames.')
return(sign_irf)
}

sign_irf <- gen_sr(loadings_full)

sign_irf[,2:ncol(sign_irf)] <- NA


TVARGIF_meboot <- function(data, ts_orig, boot = 2000){
boot = boot
store <- NULL
simdata <- NULL

for(j in 1:boot){
  simdata[[j]] <- ts_orig
}

for(i in 1:ncol(ts_orig)){
  store[[i]] <- meboot(ts_orig[,i], boot, expand.sd = T, force.clt = F, reachbnd = F, trim = 0)$ensemble
  for(j in 1:boot){
    simdata[[j]][,i] <- store[[i]][,j]
  }
}


loadings <- NULL
sign_irfs <- NULL
for(i in 1:length(simdata)){
  temp <- createfactors(simdata[[i]] %>% 
                          as_tibble() %>% 
                          mutate(date = data$date))
  loadings[[i]] <- temp$loadings
  simdata[[i]] <- temp$ts
  #Flip the signs if the loadings don't match the sign of the original loadings
  if(loadings[[i]]['RP', 'RPCGG'] > 0 & loadings[[i]]['CGG', 'RPCGG'] < 0){
    loadings[[i]]['RP', 'RPCGG'] <- -loadings[[i]]['RP', 'RPCGG']
    loadings[[i]]['CGG', 'RPCGG'] <- -loadings[[i]]['CGG', 'RPCGG']
    simdata[[i]][,'RPCGG'] <- -simdata[[i]][,'RPCGG']
  }
  sign_irfs[[i]] <- gen_sr(loadings[[i]]) 
  sign_irfs[[i]][,2:ncol(sign_irfs[[i]])] <- NA
}
rm(temp)


mods <- lapply(simdata, function(x){
  TVAR(x,
       lag = 4,
       nthresh = 1,
       mTh = 2)
  
})



out <- pbapply::pbMap(TVARGIRF_factor, tvar = mods, ts = simdata, loadings = loadings, sign_restr = sign_irfs,
                      MoreArgs = list(#sign_restr = gen_sr(loadings),
                                      #loadings = loadings_full,
                                      H = 100, K = 100)
)


for(i in 1:length(out)){
  out[[i]] <- out[[i]] %>% 
    mutate(draw = i) 
}


bind_rows(out) %>% saveRDS('data/GIRFBOOT.rds')

out2 <- NULL
for(i in 1:length(mods)){
  temp <- bind_cols(
    mods[[i]]$usedThVar %>% 
      as_tibble() %>% 
      rename(simThvar = value),
    mods[[i]]$model.specific$regime %>% 
      as_tibble() %>% 
      drop_na() %>% 
      rename(regime = value)
  ) %>% 
    mutate(gamma = mods[[i]]$model.specific$Thresh) %>% 
    mutate(draw = i) %>% 
    mutate(date = data %>% select(date) %>% slice(5:nrow(data)) %>% pull(date))
  
  out2[[i]] <- temp
}
rm(temp)

bind_rows(out2) %>% 
  group_by(draw) %>% 
  mutate(max = max(simThvar),
         x = max*ifelse(max == simThvar, 1, NA)/simThvar*regime,
         x = mean(x, na.rm = T),
         regime = ifelse(regime == x, 'Turmoil', 'Tranquil')) %>% 
  select(-max, -x) %>% 
  saveRDS('data/simThvars.rds')


}

TVARGIF_meboot(data, ts_orig)


make_thresholdsplot(tvar,data)
makegirfplot()
