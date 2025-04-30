TVARGIRF <- function(tvar, ts, H = 100, K = 500, horizon = 21, sign_restr = NULL){
  
  require(tsDyn)
  require(tidyverse)
  
  regimes <- tvar$model.specific$regime #Extract regimes from the tvar
  R <- max(regimes, na.rm = T)
  
  shock <- NULL
  for(i in 1:R){
    shock[[i]] <- tryCatch({
      #tmp <- resid(tvar)[which(regimes[!is.na(regimes)] == i), ] %>%
      tmp <- resid(tvar) %>%
        cov() %>%
        chol() %>%
        t()
      if(!is.null(sign_restr)){
        set.seed(1234)
      tmp <- sign_restr(sigma_chol = tmp,
                        sign_restr = sign_restr)
      colnames(tmp) <- colnames(ts)}
      tmp[,'R'] #/ tmp['R','R']
    }, error = function(e) {
      message(sprintf("Error in regime %d: %s. Using full residuals instead.", i, e$message))
      tmp <- resid(tvar) %>%
        cov() %>%
        chol() %>%
        t()
      if(!is.null(sign_restr)){
        set.seed(1234)
        tmp <- sign_restr(sigma_chol = tmp,
                          sign_restr = sign_restr)
        colnames(tmp) <- colnames(ts)}
      tmp[,'R'] #/ tmp['R','R']
    })
  }
  

  
  innov <- matrix(data = 0, ncol = ncol(resid(tvar)), nrow = horizon)  # Initialize matrix for the random innovations
  delta <- innov  # Matrix space for the shock
  
  
  
  #Inputs for TVAR.sim
  coef  <- tvar$coeffmat
  Thresh = tvar$model.specific$Thresh
  nthresh = tvar$model.specific$nthresh
  lag = tvar$lag
  thDelay = tvar$model.specific$thDelay
  mTh = tvar$model.specific$transCombin
  
  out <- NULL
  for(r in 1:R){
    IRF_DRAW_REGIME <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
    for(h in 1:H){
      htrack <- sample(which(regimes == r), size = 1, replace = FALSE)
      history <- ts[(htrack - lag + 1):htrack, ]
      regime <- regimes[htrack]
      delta[1, ] <- shock[[regime]]
      
      # Initialize matrix to store the averaged-out impulse responses
      IRF_DRAW <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
      colnames(IRF_DRAW) <- colnames(resid(tvar))
      
      for(k in 1:K){
        # Sample a random innovation
        restrack <- sample(which(regimes[!is.na(regimes)] == r), size = 1, replace = FALSE)
        innov[1, ] <- resid(tvar)[restrack, ] 
        
        YNULL <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
        YSHOCK <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
        
        YNULL <- TVAR.sim(B = coef,
                          Thresh = Thresh,
                          nthresh = nthresh,
                          n = horizon,
                          lag = lag,
                          thDelay = thDelay,
                          mTh = mTh,
                          starting = history,
                          innov = innov)
        
        YSHOCK <- TVAR.sim(B = coef,
                           Thresh = Thresh,
                           nthresh = nthresh,
                           n = horizon,
                           lag = lag,
                           thDelay = thDelay,
                           mTh = mTh,
                           starting = history,
                           innov = innov + delta)
        
        
        # Update the impulse response average
        IRF_DRAW <- IRF_DRAW + (YSHOCK - YNULL) / K
        
      }
      IRF_DRAW_REGIME <- IRF_DRAW_REGIME + IRF_DRAW/H
    }
    out[[r]] <- IRF_DRAW_REGIME %>% 
      as_tibble() %>% 
      rownames_to_column() %>% 
      mutate(regime = as.factor(regime)) %>% 
      rename(t = rowname) %>% 
      mutate(t = as.numeric(t)) %>% 
      gather(var, value, -t, -regime)
  }
  
  out <- bind_rows(out)
  
  regime <- tvar$model.specific$regime
  thvar <- tvar$usedThVar
  regime <- regime[which(!is.na(regime))]
  
  upper <- regime[which.max(thvar)]
  lower <- regime[which.min(thvar)]
  
  out <- out %>% 
    mutate(regime = ifelse(regime == upper, 'Turmoil', 'Tranquil'))
  
  
  out
  
}

TVARGIRF_factor <- function(tvar, ts, H = 100, K = 500, horizon = 21, sign_restr = NULL, loadings){
  
  require(tsDyn)
  require(tidyverse)
  
  regimes <- tvar$model.specific$regime #Extract regimes from the tvar
  R <- max(regimes, na.rm = T)
  
  
  shock <- NULL
  for(i in 1:R){
    shock[[i]] <- tryCatch({
      tmp <- resid(tvar)[which(regimes[!is.na(regimes)] == i), ] %>%
      #tmp <- resid(tvar) %>%
        cov() %>%
        chol() %>%
        t()
      if(!is.null(sign_restr)){
        set.seed(2000)
        tmp <- sign_restr_factor(sigma_chol = tmp,
                                 sign_restr = sign_restr,
                                 loadings = loadings)
        colnames(tmp) <- colnames(ts)}
      tmp[,'R'] #/ tmp['R','R']
    }, error = function(e) {
      message(sprintf("Error in regime %d: %s. Using full residuals instead.", i, e$message))
      tmp <- resid(tvar) %>%
        cov() %>%
        chol() %>%
        t()
      if(!is.null(sign_restr)){
        set.seed(2000)
        tmp <- sign_restr_factor(sigma_chol = tmp,
                          sign_restr = sign_restr,
                          loadings = loadings)
        colnames(tmp) <- colnames(ts)}
      tmp[,'R'] #/ tmp['R','R']
    })
  }
  
  
  
  innov <- matrix(data = 0, ncol = ncol(resid(tvar)), nrow = horizon)  # Initialize matrix for the random innovations
  delta <- innov  # Matrix space for the shock
  
  
  
  #Inputs for TVAR.sim
  coef  <- tvar$coeffmat
  Thresh = tvar$model.specific$Thresh
  nthresh = tvar$model.specific$nthresh
  lag = tvar$lag
  thDelay = tvar$model.specific$thDelay
  mTh = tvar$model.specific$transCombin
  
  out <- NULL
  for(r in 1:R){
    IRF_DRAW_REGIME <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
    for(h in 1:H){
      htrack <- sample(which(regimes == r), size = 1, replace = FALSE)
      history <- ts[(htrack - lag + 1):htrack, ]
      regime <- regimes[htrack]
      delta[1, ] <- shock[[regime]]
      
      # Initialize matrix to store the averaged-out impulse responses
      IRF_DRAW <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
      colnames(IRF_DRAW) <- colnames(resid(tvar))
      
      for(k in 1:K){
        # Sample a random innovation
        restrack <- sample(which(regimes[!is.na(regimes)] == r), size = 1, replace = FALSE)
        innov[1, ] <- resid(tvar)[restrack, ] 
        
        YNULL <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
        YSHOCK <- matrix(0, ncol = ncol(resid(tvar)), nrow = horizon)
        
        YNULL <- TVAR.sim(B = coef,
                          Thresh = Thresh,
                          nthresh = nthresh,
                          n = horizon,
                          lag = lag,
                          thDelay = thDelay,
                          mTh = mTh,
                          starting = history,
                          innov = innov)
        
        YSHOCK <- TVAR.sim(B = coef,
                           Thresh = Thresh,
                           nthresh = nthresh,
                           n = horizon,
                           lag = lag,
                           thDelay = thDelay,
                           mTh = mTh,
                           starting = history,
                           innov = innov + delta)
        
        
        # Update the impulse response average
        IRF_DRAW <- IRF_DRAW + (YSHOCK - YNULL) / K
        
      }
      IRF_DRAW_REGIME <- IRF_DRAW_REGIME + IRF_DRAW/H
    }
    out[[r]] <- IRF_DRAW_REGIME %>% 
      as_tibble() %>% 
      rownames_to_column() %>% 
      mutate(regime = as.factor(regime)) %>% 
      rename(t = rowname) %>% 
      mutate(t = as.numeric(t)) %>% 
      gather(var, value, -t, -regime)
  }
  
  out <- bind_rows(out)
  
  regime <- tvar$model.specific$regime
  thvar <- tvar$usedThVar
  regime <- regime[which(!is.na(regime))]
  
  upper <- regime[which.max(thvar)]
  lower <- regime[which.min(thvar)]
  
  out <- out %>% 
    mutate(regime = ifelse(regime == upper, 'Turmoil', 'Tranquil'))
  
  
  out
  
}


LRtests <- function(ts, mTh){
  
  require(tsDyn)
  
set.seed(2025)
test1 <- tsDyn::TVAR.LRtest(ts, 
                   lag = 4, 
                   thDelay = 4, 
                   mTh = mTh, 
                   test = '1vs',
                   nboot = 100)

set.seed(2025)
test2 <- tsDyn::TVAR.LRtest(ts, 
                   lag = 4, 
                   thDelay = 4, 
                   mTh = mTh, 
                   test = '2vs3',
                   nboot = 100)




out <- tibble(
test = c('Linear vs 1-threshold',
         'Linear vs 2-threshold',
         '1-threshold vs 2 threshold'),
LR = c(test1$LRtest.val,
       test2$LRtest.val),
P = c(test1$Pvalueboot,
      test2$Pvalueboot))

saveRDS(out, 'figures/LRtest.rds')


}

PCA_weights <- function(data){
  
  require(FactoMineR)
  require(tidyverse)
  
  PCA <- data %>% 
    as_tibble() %>% 
    select(RP, CGG) %>% 
    as.matrix() %>% 
    FactoMineR::PCA(scale.unit = F,
                    ncp = 1) 
  
  weights <- PCA$var$coord
  
  weights
}

makePCAplot <- function(data, weights) {
  require(tidyverse)
  
  pcaplot <- data %>% 
    mutate(BBIDX = weights['CGG',]*CGG + weights['RP',]*RP) %>% 
    select(date, CGG, RP, BBIDX) %>% 
    gather(var, value, -date) %>% 
    mutate(var = case_when(
      var == 'CGG'   ~ 'Credit-to-GDP gap',
      var == 'RP'    ~ 'Rent-price ratio',
      var == 'BBIDX' ~ 'PC1'
    )) %>%
    ggplot(aes(x = date, y = value, color = var, linetype = var, linewidth = var)) +
    geom_line() +
    scale_color_manual(values = c(
      'Credit-to-GDP gap' = '#4DAF4A',
      'Rent-price ratio'  = '#377EB8',
      'PC1'               = '#D73027'
    )) +
    scale_linetype_manual(values = c(
      'Credit-to-GDP gap' = 'dashed',
      'Rent-price ratio'  = 'dotted',
      'PC1'               = 'solid'
    )) +
    scale_linewidth_manual(values = c(
      'Credit-to-GDP gap' = 0.7,
      'Rent-price ratio'  = 0.7,
      'PC1'               = 1.2
    )) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "inside", 
      # Move legend slightly in from the corner
      legend.position.inside = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.background = element_rect(
        fill = "white", 
        color = "black", 
        linewidth = 0.5
      ),
      legend.title = element_blank(),
      # Make legend text smaller
      legend.text = element_text(size = 9),
      # Decrease the key size so boxes are smaller
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5))
    ) +
    labs(x = NULL, y = NULL)
  
  ggsave(
    filename = 'figures/pcaplot.pdf',
    plot = pcaplot,
    width = 6.5,
    height = 4.5
  )
}

thvar_compar <- function(ts, data, weights){
  require(tsDyn)
  require(tidyverse)
  require(fredr)
  
  fredr_set_key("cda47ae66b38ed7988c0a9c2ec80c94f")
  
  tvar1 <- tsDyn::TVAR(ts, 
                       lag = 4,
                       nthresh = 1,
                       mTh = 2)
  
  
ts2 <- ts  
ts3 <- ts
ts2[,2] = -data$RP
ts3[,2] = data$CGG

  tvar2 <- tsDyn::TVAR(ts2, 
                       lag = 4,
                       nthresh = 1,
                       mTh = 2)
  
  tvar3 <- tsDyn::TVAR(ts3, 
                       lag = 4,
                       nthresh = 1,
                       mTh = 2)
  
  
  
  regimes_data <-   data %>% 
    mutate(regime1 = tvar1$model.specific$regime,
           regime2 = tvar2$model.specific$regime,
           regime3 = tvar3$model.specific$regime) %>% 
    drop_na() %>% 
    select(date, CGG, RP, contains('regime')) %>% 
    gather(var, value, CGG, RP) %>% 
    mutate(across(contains('regime'), ~ifelse(.x == 1, 'Tranquil', 'Turmoil'))) %>% 
    mutate(var = case_when(var == 'CS' ~ 'Credit Spread',
                           var == 'HP' ~ 'Property Prices',
                           var == 'P' ~ 'Consumer Prices',
                           var == 'TS' ~ 'Term Spread',
                           var == 'HL' ~ 'Real Estate Loans',
                           var == 'MR' ~ 'Mortgage Rate',
                           var == 'R' ~ 'Fed Funds Rate',
                           var == 'Y' ~ 'Real GDP',
                           var == 'CGG' ~ 'Credit-to-GDP gap',
                           var == 'RP' ~ 'Rent-price ratio')) 
  
  
  recdum <- fredr('JHDUSRGDPBR') %>% 
    select(date, value)
  
  
  add_recession_shading <- function(p, recession_dummy) {
    recessions <- recession_dummy %>%
      mutate(recdum_lag = lag(value, default = 0),
             start = if_else(value == 1 & recdum_lag == 0, date, as.Date(NA)),
             end   = if_else(value == 0 & recdum_lag == 1, date, as.Date(NA))) %>%
      reframe(start = na.omit(start),
              end = na.omit(end)) 
    
    if (tail(recession_dummy$value, 1) == 1) {
      last_start <- recession_dummy$date[max(which(recession_dummy$value == 1))]
      recessions <- bind_rows(recessions, tibble(start = last_start, end = max(recession_dummy$date)))
    }
    
    p + 
      geom_rect(data = recessions,
                aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                fill = "gray80", alpha = 0.5,
                inherit.aes = FALSE)
  }
  
  p_pc <- regimes_data %>% 
    ggplot(aes(x = date, y = value, color = regime1, group = 1)) +
    geom_line(linewidth = .75) +
    facet_wrap(~var, scales = 'free', nrow = 1) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none", 
      # Move legend slightly in from the corner
      legend.position.inside = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.background = element_rect(
        fill = "white", 
        color = "black", 
        linewidth = 0.5
      ),
      legend.title = element_blank(),
      # Make legend text smaller
      legend.text = element_text(size = 9),
      # Decrease the key size so boxes are smaller
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5))
    ) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = c(
      'Tranquil' = '#1b9e77',
      'Turmoil' = '#c23b22'
    )) +
    labs(x = '',
         y = '',
         title = '(C) Regime indicator = PC1') +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 14))
  
  
  
  p_rp <- regimes_data %>% 
    ggplot(aes(x = date, y = value, color = regime2, group = 1)) +
    geom_line(linewidth = .75) +
    facet_wrap(~var, scales = 'free', nrow = 1) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none", 
      # Move legend slightly in from the corner
      legend.position.inside = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.background = element_rect(
        fill = "white", 
        color = "black", 
        linewidth = 0.5
      ),
      legend.title = element_blank(),
      # Make legend text smaller
      legend.text = element_text(size = 9),
      # Decrease the key size so boxes are smaller
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5))
    ) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = c(
      'Tranquil' = '#1b9e77',
      'Turmoil' = '#c23b22'
    )) +
    labs(x = '',
         y = '',
         title = '(A) Regime indicator = RP') +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 14))
  
  
  p_cgg = regimes_data %>% 
    ggplot(aes(x = date, y = value, color = regime3, group = 1)) +
    geom_line(linewidth = .75) +
    facet_wrap(~var, scales = 'free', nrow = 1) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none", 
      # Move legend slightly in from the corner
      legend.position.inside = c(0.97, 0.97),
      legend.justification = c("right", "top"),
      legend.background = element_rect(
        fill = "white", 
        color = "black", 
        linewidth = 0.5
      ),
      legend.title = element_blank(),
      # Make legend text smaller
      legend.text = element_text(size = 9),
      # Decrease the key size so boxes are smaller
      legend.key.width = unit(0.7, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5))
    ) +
    labs(x = NULL, y = NULL) +
    scale_color_manual(values = c(
      'Tranquil' = '#1b9e77',
      'Turmoil' = '#c23b22'
    )) +
    labs(x = '',
         y = '',
         title = '(B) Regime indicator = CGG') +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 14))
  
  
  
  p_pc <- add_recession_shading(p_pc, recession_dummy = recdum)
  p_rp <- add_recession_shading(p_rp, recession_dummy = recdum)
  p_cgg <- add_recession_shading(p_cgg, recession_dummy = recdum)
  
  plot <- grid.arrange(p_rp, p_cgg, p_pc, nrow = 3)
  
  
  ggsave(
    filename = 'figures/thvars.pdf',
    plot = plot,
    width = 10,
    height = 8
  )
  
  
}


draw_qi <- function(sigma_chol, sign_restr, M, i, Q, zero = zero) {
  
  if (isTRUE(zero)) {
    Q <- matrix(0, M, M)
    sel_row <- which(sign_restr[, i] == 0)
    R <- rbind(sigma_chol[sel_row, ], Q[seq_len(i - 1L), , drop = FALSE])
    qr_object <- qr(t(R))
    qr_rank <- qr_object[["rank"]]
    set <- if (qr_rank == 0) {
      seq_len(Q)
    } else {
      -seq_len(qr_rank)
    }
    N_i <- qr.Q(qr_object, complete = TRUE)[, set, drop = FALSE]
    z <- rnorm(nrow(N_i))
    N_stdn <- crossprod(N_i, z)
    q_i <- N_i %*% (N_stdn / norm(N_stdn, type = "2"))
  } else {
    if (i == 1) {
      x <- rnorm(M, 0, 1)
      q_i <- x / norm(x, type = "2")
    } else {
      x <- rnorm(M, 0, 1)
      QQ <- diag(M) - tcrossprod(Q)
      q_i <- QQ %*% x/norm(QQ %*% x, type = "2")
    }
  }
  return(q_i)
}


check_qi <- function(q_i, sigma_chol, sr_i) {
  # Compute the impulse vector from stacked matrix
  shock_vec <- sigma_chol %*% q_i
  M <- nrow(sigma_chol)  # Assuming stacked format (2M x M)
  
  sr_short <- sr_i[1:M]
  
  # Indices for non-zero sign restrictions
  sign_short_idx <- which(!is.na(sr_short) & sr_short != 0)
  
  # Zero restrictions
  zero_short_idx <- which(sr_short == 0)
  
  # Check zero restrictions: strict for short-run, tolerant for long-run
  if (any(abs(shock_vec[zero_short_idx]) > .Machine$double.eps)) return(0L)
  
  # Check sign restrictions
  svec <- sign(shock_vec)
  
  short_check <- svec[sign_short_idx]
  
  short_target <- sr_short[sign_short_idx]
  
  # All signs must either match or all be flipped
  if (identical(c(short_check), c(short_target))) {
    return(1L)
  } else if (identical(c(-short_check), c(short_target))) {
    return(-1L)
  } else {
    return(0L)
  }
}

sign_restr <- function (sigma_chol, sign_restr,
                         M=ncol(sigma_chol), sign_lim = 1000) 
{
  counter_outer <- 0L
  while (TRUE) {
    if (counter_outer > min(sign_lim^0.5)) {
      stop(paste0("No matrix fitting the sign restrictions found after ", 
                  sign_lim, " tries. Consider increasing the limit via the", 
                  "`sign_lim` or adapting the restrictions."))
    }
    counter_outer <- counter_outer + 1L
    Q <- matrix(0, M, M)
    pos_check <- logical(M)
    i <- 1L
    counter_inner <- 0L
    while (!all(pos_check)) {
      if (counter_inner > min(sign_lim^0.6)) 
        break
      counter_inner <- counter_inner + 1L
      zero <- !purrr::is_empty(which(sign_restr[,i] == 0))
      q_i <- draw_qi(sigma_chol = sigma_chol, 
                     sign_restr = sign_restr, 
                     M= M, i = i, 
                     zero = zero, Q = Q)
      sign_check <- check_qi(q_i, sigma_chol, sr_i = sign_restr[, i])
      if (sign_check != 0L) {
        pos_check[i] <- TRUE
        Q[, i] <- q_i * sign_check
        i <- i + 1L
      }
    }
    if (all(pos_check)) {
      return(sigma_chol %*% Q)
    }
  }
}

draw_qi <- function(sigma_chol, sign_restr, M, i, Q, zero = zero) {
  
  if (isTRUE(zero)) {
    Q <- matrix(0, M, M)
    sel_row <- which(sign_restr[, i] == 0)
    R <- rbind(sigma_chol[sel_row, ], Q[seq_len(i - 1L), , drop = FALSE])
    qr_object <- qr(t(R))
    qr_rank <- qr_object[["rank"]]
    set <- if (qr_rank == 0) {
      seq_len(Q)
    } else {
      -seq_len(qr_rank)
    }
    N_i <- qr.Q(qr_object, complete = TRUE)[, set, drop = FALSE]
    z <- rnorm(nrow(N_i))
    N_stdn <- crossprod(N_i, z)
    q_i <- N_i %*% (N_stdn / norm(N_stdn, type = "2"))
  } else {
    if (i == 1) {
      x <- rnorm(M, 0, 1)
      q_i <- x / norm(x, type = "2")
    } else {
      x <- rnorm(M, 0, 1)
      QQ <- diag(M) - tcrossprod(Q)
      q_i <- QQ %*% x/norm(QQ %*% x, type = "2")
    }
  }
  return(q_i)
}

check_qi_factor <- function(q_i, sigma_chol, loadings, sr_i) {
  # Compute the impulse vector from stacked matrix
  sigma_var <- loadings %*% sigma_chol
  shock_vec <- sigma_var %*% q_i
  M <- nrow(sigma_var)  # Assuming stacked format (2M x M)
  
  sr_short <- sr_i[1:M]
  
  # Indices for non-zero sign restrictions
  sign_short_idx <- which(!is.na(sr_short) & sr_short != 0)
  
  # Zero restrictions
  zero_short_idx <- which(sr_short == 0)
  
  # Check zero restrictions: strict for short-run, tolerant for long-run
  if (any(abs(shock_vec[zero_short_idx]) > .Machine$double.eps)) return(0L)
  
  # Check sign restrictions
  svec <- sign(shock_vec)
  
  short_check <- svec[sign_short_idx]
  
  short_target <- sr_short[sign_short_idx]
  
  # All signs must either match or all be flipped
  if (identical(c(short_check), c(short_target))) {
    return(1L)
  } else if (identical(c(-short_check), c(short_target))) {
    return(-1L)
  } else {
    return(0L)
  }
}

sign_restr_factor <- function (sigma_chol, sign_restr, loadings,
                               M=ncol(sigma_chol), sign_lim = 10000) 
{
  counter_outer <- 0L
  while (TRUE) {
    if (counter_outer > min(sign_lim^0.5)) {
      stop(paste0("No matrix fitting the sign restrictions found after ", 
                  sign_lim, " tries. Consider increasing the limit via the", 
                  "`sign_lim` or adapting the restrictions."))
    }
    counter_outer <- counter_outer + 1L
    Q <- matrix(0, M, M)
    pos_check <- logical(M)
    i <- 1L
    counter_inner <- 0L
    while (!all(pos_check)) {
      if (counter_inner > min(sign_lim^0.6)) 
        break
      counter_inner <- counter_inner + 1L
      zero <- !purrr::is_empty(which(sign_restr[,i] == 0))
      q_i <- draw_qi(sigma_chol = sigma_chol, 
                            sign_restr = sign_restr, 
                            M= M, i = i, 
                            zero = zero, Q = Q)
      sign_check <- check_qi_factor(q_i, sigma_chol, loadings, sr_i = sign_restr[, i])
      if (sign_check != 0L) {
        pos_check[i] <- TRUE
        Q[, i] <- q_i * sign_check
        i <- i + 1L
      }
    }
    if (all(pos_check)) {
      return(sigma_chol %*% Q)
    }
  }
}


makegirfplot <- function(){
  regime_cols <- c(
    'Tranquil' = '#1b9e77',
    'Turmoil'  = '#c23b22'
  )
  
  common_theme <- theme_minimal(base_size = 14) +
    theme(
      legend.position           = "none",
      plot.margin               = margin(t = 10, r = 10, b = 10, l = 10),
      axis.text.x               = element_text(margin = margin(t = 5)),
      axis.text.y               = element_text(margin = margin(r = 5)),
      strip.text                = element_text(size = 14, face = "bold"),
      panel.grid.minor          = element_blank(),
      plot.title                = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.title.position       = "plot"
    )
  
  # Medians:
  p_medians <-read_rds('data/GIRFBOOT.rds') %>% 
    spread(var, value) %>% 
    mutate(RP = RPCGG* loadings_full['RP', 'RPCGG'],
           CGG = RPCGG* loadings_full['CGG', 'RPCGG'],
           MR = F1*loadings_full['MR','F1'] + F2*loadings_full['MR','F2'] + F3*loadings_full['MR','F3'],
           CS = F1*loadings_full['CS','F1'] + F2*loadings_full['CS','F2'] + F3*loadings_full['CS','F3'] ,
           TS = F1*loadings_full['TS','F1'] + F2*loadings_full['TS','F2'] + F3*loadings_full['TS','F3'],
           HL = F1*loadings_full['HL','F1'] + F2*loadings_full['HL','F2'] + F3*loadings_full['HL','F3'],
           HP = F1*loadings_full['HP','F1'] + F2*loadings_full['HP','F2'] + F3*loadings_full['HP','F3'],
           P = F1*loadings_full['P','F1'] + F2*loadings_full['P','F2'] + F3*loadings_full['P','F3'],
           Y = F1*loadings_full['Y','F1'] + F2*loadings_full['Y','F2'] + F3*loadings_full['Y','F3']) %>% 
    select(-RPCGG, -F1, -F2, -F3) %>% 
    gather(var, value, -t, -regime, -draw) %>% 
    group_by(t, regime, var) %>% 
    summarize(median = median(value),
              lb = quantile(value, .16),
              ub = quantile(value, .84)) %>%
    ggplot(aes(x = t, y = median, colour = regime)) +
    geom_hline(yintercept = 0, colour = "#c23b22", linetype = "dashed") +
    geom_line(size = 0.75) +
    facet_wrap(~var, scales = "free", ncol = 2) +
    scale_colour_manual(values = regime_cols) +
    labs(title = "(a) Medians", x = NULL, y = NULL) +
    common_theme
  
  # Differences:
  p_diffs <- read_rds('data/GIRFBOOT.rds') %>% 
    spread(var, value) %>% 
    mutate(RP = RPCGG* loadings_full['RP', 'RPCGG'],
           CGG = RPCGG* loadings_full['CGG', 'RPCGG'],
           MR = F1*loadings_full['MR','F1'] + F2*loadings_full['MR','F2'] + F3*loadings_full['MR','F3'],
           CS = F1*loadings_full['CS','F1'] + F2*loadings_full['CS','F2'] + F3*loadings_full['CS','F3'] ,
           TS = F1*loadings_full['TS','F1'] + F2*loadings_full['TS','F2'] + F3*loadings_full['TS','F3'],
           HL = F1*loadings_full['HL','F1'] + F2*loadings_full['HL','F2'] + F3*loadings_full['HL','F3'],
           HP = F1*loadings_full['HP','F1'] + F2*loadings_full['HP','F2'] + F3*loadings_full['HP','F3'],
           P = F1*loadings_full['P','F1'] + F2*loadings_full['P','F2'] + F3*loadings_full['P','F3'],
           Y = F1*loadings_full['Y','F1'] + F2*loadings_full['Y','F2'] + F3*loadings_full['Y','F3']) %>% 
    select(-RPCGG, -F1, -F2, -F3) %>%
    gather(var, value, -t, -regime, -draw) %>% 
    spread(regime, value) %>% 
    mutate(diff = Turmoil - Tranquil) %>% 
    group_by(t, var) %>% 
    summarize(median = median(diff),
              lb = quantile(diff, .16),
              ub = quantile(diff, .84),
              lb2 = quantile(diff, .025),
              ub2 = quantile(diff, .975)) %>% 
    ggplot(aes(x = t, y = median)) +
    geom_ribbon(aes(ymin = lb2, ymax = ub2), fill = "grey90", alpha = .5) +
    geom_ribbon(aes(ymin = lb,  ymax = ub ), fill = "grey70", alpha = .5) +
    geom_hline(yintercept = 0, colour = "#c23b22", linetype = "dashed") +
    geom_line(size = 0.75, colour = "black") +
    facet_wrap(~var, scales = "free", ncol = 2) +
    labs(title = "(b) Differences", x = NULL, y = NULL) +
    common_theme
  
  
  girfplot <- grid.arrange(p_medians, p_diffs, ncol = 2)
  
  ggsave(
    filename = 'figures/girfplot.pdf',
    plot = girfplot,
    width = 10,
    height = 8
  )
}

make_thresholdsplot <- function(tvar, data){
  df <- read_rds('data/simThvars.rds')
  
  true_thresh <- tvar$model.specific$Thresh
  
  # custom colours for regimes
  regime_cols <- c(
    'Tranquil' = '#1b9e77',
    'Turmoil'  = '#c23b22'
  )
  
  #— Panel (a): density of γ with vertical line at true_thresh
  p_thresh <- df %>%
    distinct(draw, gamma) %>% 
    ggplot(aes(x = gamma)) +
    geom_density(fill = "grey80", colour = "black", linewidth = 0.9) +
    geom_vline(xintercept = true_thresh,
               linetype = "dashed",
               colour   = regime_cols["Turmoil"],
               size = .9) +
    annotate(
      "text",
      x      = true_thresh,
      y      = 0,
      label  = paste0("γ* = ", round(true_thresh, 4)),
      hjust  = 1.1,
      vjust  = -0.5,
      size   = 4
    ) +
    labs(
      title = "(a) Threshold density",
      x     = expression(gamma),
      y     = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid.minor   = element_blank()
    )
  
  #— Panel (b): contingency of true vs. estimated regimes
  # first compute the percentages
  tab <- tvar$model.specific$regime %>% 
    as_tibble() %>% 
    drop_na() %>% 
    rename(regime_true = value) %>% 
    mutate(date = data %>% select(date) %>% slice(5:nrow(data)) %>% pull(date)) %>% 
    inner_join(df) %>% 
    mutate (regime_true = ifelse(regime_true == 1, 'Tranquil', 'Turmoil')) %>% 
    group_by(regime_true, regime) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(pct = 100*n/sum(n)) 
  
  # then reshape and plot
  p_conf <- tab %>% 
    rename(True = regime_true, Estimate = regime) %>% 
    ggplot(aes(x = True, y = pct, fill = Estimate)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(
      aes(label = round(pct, 1)),
      position = position_dodge(width = 0.6),
      vjust = -0.5,  # controls how high above the bar the text appears
      size = 4       # controls font size
    ) +
    scale_fill_manual(values = regime_cols) +
    labs(
      title = "(b) Classification accuracy",
      x     = "True regime",
      y     = "Percent"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position  = "none"
    )
  #— Combine into a 2‐panel and print
  thresholds_plot <- grid.arrange(p_thresh, p_conf, ncol = 2)
  
  ggsave(
    filename = 'figures/thresholds.pdf',
    plot = thresholds_plot,
    width = 10,
    height = 5,
    device = cairo_pdf
  )
}


