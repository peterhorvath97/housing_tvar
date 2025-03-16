tsdata <- data %>% 
  select(-date, -RP) %>% 
  ts(frequency = 4,
     start = c(data$date %>% min %>% year,
               data$date %>% min %>% quarter))


tvpvar <- shrinkTVPVAR::shrinkTVPVAR(y = tsdata,
                                     p = 4,
                                     mod_type = 'triple',
                                     niter = 1000)


list_covmats <- function(cov, colnames){
temp2 <- NULL
for(t in 1:dim(cov)[3]){
  temp <- NULL
  for(d in 1:dim(cov)[4]){
    temp[[d]] <- cov[,,t,d]
    colnames(temp[[d]]) <- colnames
    rownames(temp[[d]]) <- colnames
  }
  temp2[[t]] <- temp
}

temp2
}
covmats_listed <- list_covmats(cov = tvpvar$Sigma,
                               colnames = colnames(tsdata))
covmats_identified <- function(covmats_listed){
T <- length(covmats_listed)
D<- length(covmats_listed[[1]])
temp2 <- NULL
for(t in 1:T){
  temp <- NULL
  for(d in 1:D){
    covmat <- covmats_listed[[t]][[d]]
    covmat <- chol(covmat)
    covmat <- t(covmat)
    covmat <- covmat/diag(covmat)
    temp[[d]] <- covmat 
  }
  temp2[[t]] <- temp
}
temp2 
}
cholesky_listed <- covmats_identified(covmats_listed)

list_coefs <- function(coef, colnames){
L <- dim(coef)[3]  
T <- dim(coef)[4]
D <- dim(coef)[5]

temp3 <- NULL
for(t in 1:T){
  temp2 <- NULL
  for(d in 1:D){
    temp <- NULL
    for(l in 1:L){
     temp[[l]] <- coef[,,l,t,d]
     rownames(temp[[l]]) <- colnames
     colnames(temp[[l]]) <- paste0(colnames, ".l", l)
    }
    temp2[[d]] <- temp
  }  
  temp3[[t]] <- temp2
}
temp3
}
coefs_listed <- list_coefs(coef = tvpvar$beta,
                           colnames = colnames(tsdata))




compute_irf <- function(var_coeffs, struct_shocks, h = 20) {
  # Extract variable names from structural shocks
  var_names <- colnames(struct_shocks)
  
  # Get the number of variables
  num_vars <- ncol(struct_shocks)
  num_lags <- length(var_coeffs)
  
  # Initialize an array to store impulse responses (horizon x variables x shocks)
  irf <- array(0, dim = c(h, num_vars, num_vars), 
               dimnames = list(1:h, var_names, paste0("", var_names)))
  
  # Loop over each variable to compute its response to each shock
  for (shock_var in 1:num_vars) {
    # Initialize a response matrix
    response <- matrix(0, nrow = h, ncol = num_vars)
    
    # Initial impact is just the shock itself
    response[1, ] <- struct_shocks[,shock_var] 
    
    # Apply VAR dynamics
    for (t in 2:h) {
      for (lag in 1:min(num_lags, t-1)) {
        response[t, ] <- response[t, ] + var_coeffs[[lag]] %*% response[t - lag, ]
      }
    }
    
    # Store the impulse response for the given shock
    irf[, , shock_var] <- response
  }
  
  return(irf)
}
sim_irfs <- function(coefs_listed, cholesky_listed){
T <- length(cholesky_listed)
D <- length(cholesky_listed[[1]])

temp2 <- NULL
for(t in 1:T){
  temp <- NULL
  for(d in 1:D){
    temp[[d]] <- compute_irf(var_coeffs = coefs_listed[[t]][[d]],
                             struct_shocks = cholesky_listed[[t]][[d]],
                             h = 21) %>% 
      as.data.frame.table() %>% 
      as_tibble() %>% 
      rename(horizon = Var1,
             var = Var2,
             shock = Var3,
             value = Freq) %>% 
      mutate(horizon = as.numeric(horizon)-1,
             time = t,
             draw = d)  
    
  }
  temp <- bind_rows(temp)
  temp2[[t]] <- temp
}
temp2 <- bind_rows(temp2)
temp2
}
irfs <- sim_irfs(coefs_listed = coefs_listed,
                 cholesky_listed = cholesky_listed)



irfs_fin <- irfs %>% 
  group_by(time, draw) %>% 
  mutate(is_unstable = any(abs(value) > 10)) %>% 
  ungroup %>% 
  filter(!is_unstable) %>% 
  select(-is_unstable) %>% 
  group_by_at(vars(-value, -draw)) %>% 
  summarize(irf = median(value),
            ub1 = quantile(value, probs = 0.84),
            ub2 = quantile(value, probs = 0.975),
            lb1 = quantile(value, probs = 0.16),
            lb2 = quantile(value, probs = 0.025)) %>% 
  ungroup()


irfs_fin <- data %>% 
  slice(5:nrow(data)) %>% 
  select(date) %>% 
  rownames_to_column() %>% 
  mutate(rowname = as.numeric(rowname)) %>% 
  inner_join(irfs_fin, by = c('rowname' = 'time')) %>% 
  select(-rowname, time = date)




irfs_fin %>% 
  filter(shock == 'R') %>%
  filter(time == '2015-10-01') %>% 
  mutate(time = as.factor(time)) %>% 
  ggplot(aes(x = horizon, y = irf, group = time)) +
  geom_line(aes(color = time), linewidth = .8) +
  geom_ribbon(aes(ymin = lb1, ymax = ub1), alpha = .2, fill = 'blue') +
  #geom_ribbon(aes(ymin = lb2, ymax = ub2), alpha = .1, fill = 'blue') +
  geom_hline(aes(yintercept = 0), color = 'red') +
  facet_wrap(~var, ncol = 2, scales = 'free') +
  theme_minimal() +
  labs(x = '', y = '')


