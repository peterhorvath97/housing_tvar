data_ts <- lapply(data, function(data){data %>% select(-date) %>% ts()})


defineBmat <- function(varest){
  n <- varest %>% 
    resid() %>% 
    ncol()
  
  mat <- matrix(nrow = n, ncol = n, 0)
  rownames(mat) <- colnames(varest %>% resid())
  colnames(mat) <- colnames(varest %>% resid())
  mat[lower.tri(mat)] <- NA
  mat[,'R'] <- 0
  mat[,'M30'] <- 0
  mat['M30','R'] <- NA
  mat['HLOAN','R'] <- NA
  mat['HLOAN','M30'] <- NA
  #mat['R', 'Y'] <- NA
  #mat['R', 'P'] <- NA
  if(max(str_detect(colnames(mat), '\\bC\\b')) == 1){
    mat[,'C'] <- 0
  }
  if(max(str_detect(colnames(mat), '\\bInv\\b')) == 1){
    mat[,'Inv'] <- 0
  }
  diag(mat) <- 1
  mat
}
compute_impulse_responses <- function(var_model, shocks, horizon = 10) {
  # Extract coefficient matrices from the VAR model
  k <- var_model$K  # Number of endogenous variables
  p <- var_model$p  # Number of lags
  
  # Flatten coefficient matrix to a list of lag matrices
  coef_matrices <- array(0, dim = c(k, k, p))
  for (i in 1:p) {
    coef_matrices[, , i] <- vars::Acoef(var_model)[[i]]
  }
  
  # Initialize impulse response matrix
  IR <- matrix(0, nrow = p+horizon, ncol = k)
  
  # Set the initial impulse (impact period)
  IR[p+1, ] <- shocks
  
  
  # Compute impulse responses recursively
  for (t in (p+2):(p+horizon)) {
    response <- matrix(0, nrow = k, ncol = 1)
    for (j in 1:p) {
      response <- response + coef_matrices[, , j] %*% IR[t-j, ]
    }
    IR[t, ] <- response
  }
  
  IR <- IR[(p+1):(horizon+p), ]
  colnames(IR) <- names(shocks)
  IR <- as_tibble(IR) %>% 
    rownames_to_column() %>% 
    rename(t = rowname) %>% 
    mutate(t = as.numeric(t)) %>% 
    gather(var, value, -t)
  IR <- IR %>% 
    mutate(var = factor(var,levels = unique(IR$var)))
  
  IR
  
}
est.lin <- function(estdata, p){
  
  require(tidyverse)
  
  VAR <- estdata %>% 
    vars::VAR(p = p)
  
  SVAR <- vars::SVAR(VAR,
                     Bmat = defineBmat(VAR),
                     estmethod = 'scoring')
  
  IRF <- vars::irf(SVAR, 
                   impulse = 'R',
                   response = colnames(SVAR$B),
                   n.ahead = 5*p+1,
                   ci = .68) 
  
  
  irf <- IRF$irf$R %>% 
    as_tibble() %>% 
    mutate(t = row_number()) %>% 
    gather(var, irf, -t) 
  
  ub <- IRF$Upper$R %>% 
    as_tibble() %>% 
    mutate(t = row_number()) %>% 
    gather(var, ub, -t) 
  
  lb <- IRF$Lower$R %>% 
    as_tibble() %>% 
    mutate(t = row_number()) %>% 
    gather(var, lb, -t) 
  
  IRF <- inner_join(irf, ub) %>% 
    inner_join(lb)
  
  
  IRF <- IRF %>% 
    mutate(var = factor(var,levels = unique(IRF$var)))
  
  IRF2 <- compute_impulse_responses(var_model = VAR, shocks = SVAR$B[,'R'], 5*p+1)
  
  list(VAR = VAR, SVAR = SVAR, IRF = IRF, IRF2 = IRF2)
}

est_lin <- list(Q = est.lin(data_ts$dataQ, 4),
                Qplus = est.lin(data_ts$dataQplus, 4),
                M = est.lin(data_ts$dataM, 12))

est_lin$M$IRF %>% 
  ggplot(aes(x = t, y = irf)) +
  geom_line() +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .2) +
  geom_hline(linetype = 'dashed', color = 'red', yintercept = 0) +
  facet_wrap(~var, scales = 'free') +
  theme_minimal() +
  labs(x = '', 
       y = '') +
  geom_line(data = est_lin$M$IRF2, aes(x = t, y = value), color = 'red', linewidth = 2)


Acoef_tvar <- function(var_model, regime, p){
  # Extract matrix from your object (assuming it's already stored in Bdown)
  coefmat <- as.matrix(var_model$coefficients[[regime]])
  
  # Get the number of lags from column names
  col_names <- colnames(coefmat)
  lags = p
  
  
  # Number of equations
  num_equations <- nrow(coefmat)
  
  # Initialize an empty list to store the lag matrices
  lag_matrices <- vector("list", lags)
  
  # Loop through lags and extract sub-matrices
  for (lag in 1:lags) {
    lag_cols <- grep(paste0("-", lag, "$"), col_names)
    lag_matrices[[lag]] <- coefmat[, lag_cols, drop = FALSE]
  }
  
  # Optionally, convert to a 3D array if needed
  lag_array <- array(unlist(lag_matrices), dim = c(num_equations, ncol(lag_matrices[[1]]), lags))
  
  # Display the result
  lag_array
}
compute_impulse_responses_tvar <- function(var_model, shocks, horizon = 10) {
  # Extract coefficient matrices from the VAR model
  r <- length(var_model$nobs_regimes)
  k <- var_model$k  # Number of endogenous variables
  p <- var_model$lag  # Number of lags
  
  IRS <- NULL
  for(regime in 1:r){
    # Flatten coefficient matrix to a list of lag matrices
    coef_matrices <- array(0, dim = c(k, k, p))
    for (i in 1:p) {
      coef_matrices[, , i] <- Acoef_tvar(var_model, regime = regime, p = p)[,,i]
    }
    
    # Initialize impulse response matrix
    IR <- matrix(0, nrow = p+horizon, ncol = k)
    
    # Set the initial impulse (impact period)
    IR[p+1, ] <- shocks
    
    
    # Compute impulse responses recursively
    for (t in (p+2):(p+horizon)) {
      response <- matrix(0, nrow = k, ncol = 1)
      for (j in 1:p) {
        response <- response + coef_matrices[, , j] %*% IR[t-j, ]
      }
      IR[t, ] <- response
    }
    
    IR <- IR[(p+1):(horizon+p), ]
    colnames(IR) <- names(shocks)
    IR <- as_tibble(IR) %>% 
      rownames_to_column() %>% 
      rename(t = rowname) %>% 
      mutate(t = as.numeric(t)) %>% 
      gather(var, value, -t)
    IR <- IR %>% 
      mutate(var = factor(var,levels = unique(IR$var)),
             regime = as_factor(regime))
    
    IRS[[regime]] <- IR
    
  }
  
  IRS <- bind_rows(IRS) 
  
  IRS
  
}
GIRF <- function(tvar, horizon, reps, hists, shock, seed, data = data,
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
  
  for(k in start:tvar$T){
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
        Y_null <- tsDyn::TVAR.sim(B = tvar$coeffmat,
                                  Thresh = tvar$model.specific$Thresh,
                                  nthresh = tvar$model.specific$nthresh,
                                  n = horizon,
                                  lag = tvar$lag,
                                  thDelay = tvar$model.specific$thDelay,
                                  mTh = 5,
                                  starting = data[((k-tvar$lag+1):k),],
                                  innov = resmat)
        
        Y_shock <- tsDyn::TVAR.sim(B = tvar$coeffmat,
                                   Thresh = tvar$model.specific$Thresh,
                                   nthresh = tvar$model.specific$nthresh,
                                   n = horizon,
                                   lag = tvar$lag,
                                   thDelay = tvar$model.specific$thDelay,
                                   mTh = 5,
                                   starting = data[((k-tvar$lag+1):k),],
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



est.tvar <- function(estdata, p){
mTh = which(str_detect(colnames(estdata), 'HPRICE') == T)
#tsDyn::TVAR.LRtest(data_ts$dataM, lag = p, mTh = mTh, nboot = 100, plot = T)

T2 <- tsDyn::TVAR(estdata, 
                  lag = p, 
                  nthresh = 1, 
                  mTh = mTh, 
                  model = "TAR", 
                  max.iter = 100) 
T3 <- tsDyn::TVAR(estdata, 
                  lag = p, 
                  nthresh = 2, 
                  mTh = mTh, 
                  model = "TAR", 
                  max.iter = 100) 



IND2 <- T2$model.specific$regime
IND3 <- T3$model.specific$regime

source('codes/SVAR2.R')

ST2 <- SVAR2(T2, estdata,
             Bmat = defineBmat(T2),
             estmethod = 'scoring')

ST3 <- SVAR2(T3, estdata,
             Bmat = defineBmat(T3),
             estmethod = 'scoring')

IRFS_COEF2 <- compute_impulse_responses_tvar(var_model = T2, shocks = ST2$B[,'R'], 5*p+1)

IRFS_COEF3 <- compute_impulse_responses_tvar(var_model = T3, shocks = ST3$B[,'R'], 5*p+1)



GIRF2_UNGROUP <- GIRF(tvar = T2,
                      data = estdata,
                      horizon = 5*p+1,
                      reps = 10,
                      hists = 10,
                      shock = t(ST2$B[,'R']),
                      seed = 1234,
                      aggregate = "median",
                      prob = 0.16)

GIRF3_UNGROUP <- GIRF(tvar = T3,
                      data = estdata,
                      horizon = 5*p+1,
                      reps = 10,
                      hists = 10,
                      shock = t(ST3$B[,'R']),
                      seed = 1234,
                      aggregate = "median",
                      prob = 0.16)

GIRF2_GROUP <- GIRF(tvar = T2,
                    data = estdata,
                    horizon = 5*p+1,
                    reps = 10,
                    hists = 10,
                    shock = t(ST2$B[,'R']),
                    seed = 1234,
                    aggregate = "median",
                    prob = 0.16,
                    group = TRUE)

GIRF3_GROUP <- GIRF(tvar = T3,
                    data = estdata,
                    horizon = 5*p+1,
                    reps = 10,
                    hists = 10,
                    shock = t(ST3$B[,'R']),
                    seed = 1234,
                    aggregate = "median",
                    prob = 0.16,
                    group = TRUE)

list(R2 = list(TVAR = T2,
               indicator = IND2,
               STVAR = ST2,
               IRF_COEF = IRFS_COEF2,
               GIRF_UNGROUP = GIRF2_UNGROUP,
               GIRF_GROUP = GIRF2_GROUP),
     R3 = list(TVAR = T3,
               indicator = IND3,
               STVAR = ST3,
               IRF_COEF = IRFS_COEF3,
               GIRF_UNGROUP = GIRF3_UNGROUP,
               GIRF_GROUP = GIRF3_GROUP)
)
}
est_tvar <- list(Q = est.tvar(data_ts$dataQ, 4),
                Qplus = est.tvar(data_ts$dataQplus, 4),
                M = est.tvar(data_ts$dataM, 12))



data$dataM %>% 
  bind_cols(T3$model.specific$regime) %>% 
  rename(indicator = ...8) %>% 
  drop_na() %>% 
  mutate(indicator = as.factor(indicator)) %>% 
  gather(var, value, -date, -indicator) %>% 
  ggplot(aes(x = date, y = value, color = indicator, group = 1)) +
  geom_line(linewidth = 1) +
  facet_wrap(~var, scales = "free") + 
  theme_minimal() 



#LR 2 vs 3
1 - pchisq(2*(logLik(T3)-logLik(T2)), df = T3$k*T3$lag*length(T3$nobs_regimes)-T2$k*T2$lag*length(T2$nobs_regimes))


est_tvar$M$R3$GIRF_GROUP %>% 
  mutate(regime = as.factor(regime)) %>% 
  ggplot(aes(x = t, y = response, color = regime)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  geom_hline(linetype = 'dashed', color = 'red', yintercept = 0) +
  facet_wrap(regime~variable, scales = 'free') +
  theme_minimal() +
  labs(x = '', 
       y = '') 

