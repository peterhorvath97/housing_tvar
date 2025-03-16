library(bvars)

tsdata <- data %>% 
  select(-date) %>% 
  ts(frequency = 4,
     start = c(data$date %>% year %>% min,
               data$date %>% quarter %>% min))

mn <- set_prior_minnesota(tsdata, nolags = 4)
ui <- set_prior_uninformative(tsdata, nolags = 4)

bv1 <- bvar(tsdata, 
            mn, 
            nreps = 1000, 
            burnin = 100)
bv2 <- bvar(tsdata, 
            ui, 
            nreps = 1000, 
            burnin = 100)



btv1 <- tvar(tsdata, 
             mn, 
             nreps = 1000, 
             burnin = 100, 
             thVar = which(colnames(tsdata) == 'RP'), 
             stabletest = F)
btv2 <- tvar(tsdata, 
             ui, 
             nreps = 1000, 
             burnin = 100, 
             thVar = which(colnames(tsdata) == 'RP'), 
             stabletest = F)



bv1_irf <- bvars::irf(bv1, 
                      id_obj = set_identification_cholesky(), 
                      nhor = 20, 
                      ncores = parallel::detectCores()-1, 
                      irfquantiles = c(0.05, 0.95, 0.32, 0.68))

bv2_irf <- bvars::irf(bv2, 
                      id_obj = set_identification_cholesky(), 
                      nhor = 20, 
                      ncores = parallel::detectCores()-1, 
                      irfquantiles = c(0.05, 0.95, 0.32, 0.68))





btv1_irf <- bvars::irf(btv1, 
                       id_obj = set_identification_cholesky(), 
                       nhor = 20, 
                       ncores = parallel::detectCores()-1, 
                       irfquantiles = c(0.05, 0.95, 0.32, 0.68))

btv2_irf <- bvars::irf(btv2, 
                       id_obj = set_identification_cholesky(), 
                       nhor = 20, 
                       ncores = parallel::detectCores()-1, 
                       irfquantiles = c(0.05, 0.95, 0.32, 0.68))
