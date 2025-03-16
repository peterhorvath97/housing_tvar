GIRF2 <- function(tvar, shock, horizon = 20, H = 200, R = 500, restrict.to = NA) {
  if (length(shock) != tvar$k) {
    stop(paste0("Your shock vector has the wrong length. Should be length ", tvar$k, " (the number of variables in your TVAR), but you passed in a ", class(shock), " with length ", length(shock)))
  }
  
  if (length(tvar$model.specific$transCombin) == 0) {
    stop("Simulations of models with external transition variables (argument thVar in TVAR) not supported",
         call. = FALSE)
  }
  
  data <- tvar$model[, 1:tvar$k]
  
  # Split the residuals by regime
  resdis <- list()
  for (r in 1:tvar$model.specific$nreg) {
    resdis[[r]] <- tvar$residuals[r == tvar$model.specific$regime[(1+tvar$lag):length(tvar$model.specific$regime)], ]
  }
  
  Y <- matrix(0, nrow = horizon, ncol = tvar$k)
  
  pb <- progress::progress_bar$new(total = H)
  
  for (h in 1:H) {
    pb$tick()
    
    # Find a history
    got.regime <- FALSE
    while (!got.regime) {
      start <- sample(tvar$t + 1,1)
      history <- data[start:(start+tvar$lag-1), ] # lower case omega t-1 in KPP notation
      r <- tvar$model.specific$regime[start + tvar$lag]
      if (is.na(restrict.to) || r == restrict.to) {
        got.regime <- TRUE
      }
    }
    
    # Now repeatedly simulate with and without shock
    Y_shock <- matrix(0, nrow = horizon, ncol = tvar$k)
    Y_base <- matrix(0, nrow = horizon, ncol = tvar$k)
    for (i in 1:R) {
      Y_shock <- Y_shock + GIRF.sim(tvar, history, horizon, shock, resdis)
      Y_base <- Y_base + GIRF.sim(tvar, history, horizon, NULL, resdis)
    }
    
    # Add the results from the history to the accumulator
    Y <- Y + (1/R)*(Y_shock - Y_base)
  }
  
  # Scale by number of histories
  Y <- (1/H)*Y
  colnames(Y) <- colnames(tvar$model)[1:tvar$k]
  Y <- tibble::as_tibble(Y)
  return(structure(list(
    responses = Y,
    H = H,
    R = R,
    shock = shock,
    tvar_name = deparse(substitute(tvar))
  ),
  class = "tvarGIRF"))
}

getregime <- function(tvar, input) {
  threshval <- 0
  for (i in 1:ncol(input)) {
    threshval <- threshval + input[1, i]*tvar$model.specific$transCombin[i]
  }
  r <- 1 + sum(threshval > tvar$model.specific$Thres)
  return(r)
}

GIRF.sim <- function(tvar, history, horizon, shock, resdis) {
  r <- tvar$model.specific$regime[start + tvar$lag]
  if(is.null(shock)) {
    # no imposed shock, so bootstrap one
    s <- sample(nrow(resdis[[r]]), size=1)
    shock <- resdis[[r]][s, ]
  }
  
  Y <- matrix(0, nrow = horizon, ncol = tvar$k)
  Y[1, ] <- sim.advance(tvar, history, shock)
  if (nrow(history) > 1) {
    history <- rbind(history[1:(nrow(history)-1), , drop = FALSE], Y[1, ])
  } else {
    history <- matrix(Y[1, ], nrow = 1, ncol = tvar$k)
  }
  shocklist <- shock
  for (t in 2:horizon) {
    r <- tvar$model.specific$regime[start + tvar$lag]
    # Sample a new shock
    s <- sample(nrow(resdis[[r]]), size=1)
    shock <- resdis[[r]][s, ]
    
    shocklist <- rbind(shocklist, shock)
    
    Y[t, ] <- sim.advance(tvar, history, shock)
    if (nrow(history) > 1) {
      history <- rbind(history[2:nrow(history), , drop = FALSE], Y[t, ])
    } else {
      history <- matrix(Y[t, ], nrow = 1, ncol = tvar$k)
    }
  }
  return(Y)
}

sim.advance <- function(tvar, history, v) {
  return(tsDyn::TVAR.sim(B = tvar$coeffmat,
                         Thresh = tvar$model.specific$Thresh,
                         nthres = tvar$model.specific$nthresh,
                         n = 1, lag = tvar$lag, include = tvar$include,
                         thDelay = tvar$model.specific$thDelay, mTh = which(tvar$model.specific$transCombin == 1),
                         starting = history, innov = matrix(data = v, nrow = 1, ncol = tvar$k)))
}

