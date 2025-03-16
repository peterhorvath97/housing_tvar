
plot.bvirf <- function(x,...){
  
  # declare variable
  irfObj <- x
  nLength <- irfObj$no_variables
  pltList <- list()
  
  if(is.null(irfObj$varnames)){
    
    varnames <- as.character(seq(1:nLength))
    
  }else{varnames <- irfObj$varnames}
  
  for(ii in 1:nLength){
    for(jj in 1:nLength){
      irf1     <- irfObj$irfdraws[ii,jj,,1]
      irfUpper <- irfObj$irfdraws[ii,jj,,2]
      irfLower <- irfObj$irfdraws[ii,jj,,3]
      irfLength <- length(irf1)
      
      # Put all the information into a data frame
      irfDf <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper, Lower = irfLower) %>% 
        as_tibble() %>% 
        mutate(shock = varnames[ii],
               var = varnames[jj]) 
      
      # Store all plots in a list
      pltList[[(jj-1)*nLength+ii]] <- irfDf
    }
  }
  
bind_rows(pltList) %>% 
  rename(t = x) 

}



plot.msirf <- function(x,...){
  
  irfObj <- x
  nLength <- irfObj$no_variables
  outlist <- list()
  
  if(is.null(irfObj$varnames)){
    
    varnames <- as.character(seq(1:nLength))
    
  }else{varnames <- irfObj$varnames}
  
  # Plot over all regimes
  
  for(kk in 1:irfObj$noregimes){
    pltList <- list()
    # Plot over all regimes
    for(ii in 1:nLength){
      for(jj in 1:nLength){
        irf1     <- irfObj$irfdraws[ii,jj,,1,kk]
        irfUpper <- irfObj$irfdraws[ii,jj,,2,kk]
        irfLower <- irfObj$irfdraws[ii,jj,,3,kk]
        irfLength <- length(irf1)
        
        # Put all the information into a data frame
        irfDf <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper, Lower = irfLower) %>% 
          as_tibble() %>% 
          mutate(shock = varnames[ii],
                 var = varnames[jj],
                 regime = as.factor(kk)) 

        pltList[[(jj-1)*nLength+ii]] <- irfDf
        
      }
    }
outlist[[kk]] <- bind_rows(pltList)
    
  }
  
  bind_rows(outlist) %>% 
    rename(t = x)
  
  
}



plot.tvirf <- function(x,...){
  
  # Initialize list to store impulse-response functions
  nLength <- dim(x$irf)[1]
  pltListReg1 <- list()
  pltListReg2 <- list()
  
  if(is.null(x$varnames)){
    
    varnames <- as.character(seq(1:nLength))
    
  }else{ varnames <- x$varnames}
  
  for(ii in 1:nLength){
    for(jj in 1:nLength){
      
      irf1 <- x$irf[ii,jj,,1,1]
      irf2 <- x$irf[ii,jj,,2,1]
      
      irfLower1 <- x$irf[ii,jj,,1,2]
      irfLower2 <- x$irf[ii,jj,,2,2]
      
      irfUpper1 <- x$irf[ii,jj,,1,3]
      irfUpper2 <- x$irf[ii,jj,,2,3]
      
      irfLength <- length(irf1)
      
      irfDf1 <- data.frame(x = seq(1:irfLength), irf = irf1, Upper = irfUpper1, Lower = irfLower1) %>% 
        as_tibble() %>% 
        mutate(shock = varnames[ii],
               var = varnames[jj],
               regime = 1) 
      irfDf2 <- data.frame(x = seq(1:irfLength), irf = irf2, Upper = irfUpper2, Lower = irfLower2) %>% 
        as_tibble() %>% 
        mutate(shock = varnames[ii],
               var = varnames[jj],
               regime = 2) 
      
   
      pltListReg1[[(jj-1)*nLength+ii]] <- irfDf1
      pltListReg2[[(jj-1)*nLength+ii]] <- irfDf2
      
    }
  }
  bind_rows(
    bind_rows(pltListReg1),
    bind_rows(pltListReg2)
  ) %>% 
    rename(t = x) %>% 
    mutate(regime = as.factor(regime)) %>% 
    filter(shock == 2) %>%
    ggplot(aes(x = t, y = irf, ymin = Lower, ymax = Upper, color = regime,
               fill = regime, group = regime)) +
    geom_line() + 
    geom_ribbon(alpha = .2) +
    geom_hline(aes(yintercept = 0), color = 'red') +
    facet_wrap(~var, scales = 'free') +
    theme_minimal()
  
}

forecast.tvar <- function(obj,forecastHorizon,interval= c(0.05,0.95),...){
  
  
  # Preliminary Calculations
  nVariables       <- dim(obj$data_info$data)[2] # Number of variables
  nLength          <- dim(obj$data_info$data)[1] # Length of time series
  nLags            <- obj$general_info$nolags  # Number of lags
  nForecasts       <- dim(obj$mcmc_draws$Alpha)[4] # Number of forecasts, depends on the number sampled posteriors
  
  nl1 <- max(obj$general_info$thMax,obj$general_info$nolags)
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nl1,nVariables,nForecasts)) # Matrix to storage forecasts
  
  print(dim(mForecastStorage))
  
  
  if(is.ts(obj$data_info$data)){
    
    tsStart          <- start(obj$data_info$data)
    tsFrequency      <- frequency(obj$data_info$data)
    
  }
  
  for(ii in 1:nForecasts){
    
    mForecastStorage[1:nl1,,ii] <- (obj$data_info$data[nLength:(nLength-nl1+1),])
    
    for(jj in 1:forecastHorizon){
      
      nStart <- jj
      nEnd   <- jj+nLags-1
      y <- mForecastStorage[nStart:nEnd,,ii]
      thCheck <- mForecastStorage[jj + nl1 - obj$mcmc_draws$deldraws[ii] ,obj$general_info$thVar, ii]
      
      if(thCheck < obj$mcmc_draws$tardraws[ii]){
        
        # First Regime
        Alpha <- obj$mcmc_draws$Alpha[,,1,ii]
        Sigma <- obj$mcmc_draws$Sigma[,,1,ii]
        
      }else{
        
        # Second Regime
        Alpha <- obj$mcmc_draws$Alpha[,,2,ii]
        Sigma <- obj$mcmc_draws$Sigma[,,2,ii]
        
      }
      
      randDraw <- stats::rnorm(nVariables) %*% t(chol(Sigma))
      
      if(obj$general_info$intercept){
        
        y <- c(1,t(y))
        tempForecast <- y %*% Alpha + randDraw
        
      }else{
        y <- t(y)
        tempForecast <- y %*% Alpha + randDraw
      }
      
      # Storing the forecast
      mForecastStorage[jj+nl1,,ii] <- tempForecast
      
    }
    
  }
  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nl1),,]
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))
  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){
      
      forecastMean[jj,ii]  <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=max(interval))
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=min(interval))
      
    }
  }
  
  #forecastFinalMean  <- rbind(as.matrix(tvarObj$mydata),forecastMean)
  #forecastFinalUpper <- rbind(as.matrix(tvarObj$mydata),forecastUpper)
  #forecastFinalLower <- rbind(as.matrix(tvarObj$mydata),forecastLower)
  OriginalPath       <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean  <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalUpper <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalLower <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  
  OriginalPath[1:nLength,] <-obj$data_info$data
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),] <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower
  
  if(is.ts(obj$data_info$data)){
    
    forecastFinalMean  <- ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- ts(OriginalPath,start=tsStart,frequency=tsFrequency)
    
  }
  
  colnames(forecastFinalMean)  <- colnames(obj$data_info$data)
  colnames(forecastFinalUpper) <- colnames(obj$data_info$data)
  colnames(forecastFinalLower) <- colnames(obj$data_info$data)
  colnames(OriginalPath)       <- colnames(obj$data_info$data)
  
  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fcbvar")
  return(retList)
  
  
}