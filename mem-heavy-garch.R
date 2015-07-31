install.packages(с("dplyr", "xts"))
require(dplyr)
require(xts)

##---------------------------------------------------------------------------
## Import HEAVY estimators from EViews through CSV
##---------------------------------------------------------------------------

## heavyStaticEstim      <- list()
## heavyStaticEstim$sp   <- read.csv("./coef_2008/heavy_sp_coef.csv")
#heavyStaticEstim$stox <- read.csv("./coef_2008/heavy_stox_coef.csv")

heavy.static.coef.sp <- read.csv("./coef_2008/heavy_coef_sp.csv")
heavy.static.coef.rs <- read.csv("./coef_2008/heavy_coef_rs.csv")

##--------------------------------------------------------------------------
## Get HEAVY coefficients from heavyStaticEstim
##--------------------------------------------------------------------------


## joinHeavyEstim  <- join_all(heavyStaticEstim, by = "Variable")
## coefStaticHeavy <- joinHeavyEstim[ ,colnames(joinHeavyEstim) == "Coef"]
## rownames(coefStaticHeavy) <- joinHeavyEstim$Variable
## colnames(coefStaticHeavy) <- names(heavyStaticEstim)


##
## TODO: add discription
##

heavy.static.test.sp <-
    cbind(dtaTest[ ,"SPX2.openprice"],
          dtaTest[ ,"SPX2.closeprice"],
          dtaTest[ ,"SPX2.rk"]) %>%
    na.omit() %>%
    data.frame() %>%
    mutate(rcto  = log(SPX2.closeprice / SPX2.openprice),
           garch = VarStaticForc(init.var = dtaEViewsSP$rcto %>% var(),
                                 const    = heavy.static.coef.sp[1,  "Coef"],
                                 coef     = heavy.static.coef.sp[3:4,"Coef"],
                                 coef.var = heavy.static.coef.sp[2,  "Coef"],
                                 data     = cbind(SPX2.rk, rcto) %>% head(-1)),
           mse   = (SPX2.rk - garch)^2,
           qlike = log(garch) - SPX2.rk / garch)


heavy.static.test.rs <-
    cbind(dtaTest[ ,"RUT2.openprice"],
          dtaTest[ ,"RUT2.closeprice"],
          dtaTest[ ,"RUT2.rk"]) %>%
    na.omit() %>%
    data.frame() %>%
    mutate(rcto   = log(RUT2.closeprice / RUT2.openprice),
           ind.rk = ifelse(rcto < 0, 1, 0) * RUT2.rk,
           garch  = VarStaticForc(init.var = dtaEViewsRUS$rcto %>% var(),
                                  const    = heavy.static.coef.rs[1,  "Coef"],
                                  coef     = heavy.static.coef.rs[3:4,"Coef"],
                                  coef.var = heavy.static.coef.rs[2,  "Coef"],
                                  data     = cbind(ind.rk, rcto) %>% head(-1)),
           mse   = (RUT2.rk - garch)^2,
           qlike = log(garch) - RUT2.rk / garch)


VarStaticForc <- function(init.var, const, coef, coef.var, data){

# A recursive computation of one-step-ahead static forecasts.

# @operator %*%    -- matrix multiplication
# @param init.var  -- initial value of variance
# @param const     -- estimator of GARCH constant
# @param coef      -- estimators of GARCH coefficients
# @param coef.var  -- estimator of conditional variance
# @param data      -- data
# @return          -- vector of one-step-ahead static forecasts

    if (nrow(data) == 0){
        init.var
    } else {
        tmp <- VarStaticForc(init.var, const, coef, coef.var, head(data, -1))
        c(tmp, const + tail(data, 1) %*% coef + coef.var * tail(tmp, 1))
  }
}

##------------------------------------------------------------------------
## Check static one-step-ahead forecast
## ----------------------------------------------------------------------


CheckStaticForc <- function(forc, prev.var, const, coef, coef.var, prev.data){

   forc == const + as.numeric(prev.data)%*%coef + coef.var*prev.var
}


CheckStaticForc(heavy.static.test.sp[4,5],
                heavy.static.test.sp[3,5],
                heavy.static.coef.sp[1,"Coef"],
                heavy.static.coef.sp[3:4,"Coef"],
                heavy.static.coef.sp[2,"Coef"],
                heavy.static.test.sp[3,3:4])

CheckStaticForc(heavy.static.test.rs[6,6],
                heavy.static.test.rs[5,6],
                heavy.static.coef.rs[1,"Coef"],
                heavy.static.coef.rs[3:4,"Coef"],
                heavy.static.coef.rs[2,"Coef"],
                c(heavy.static.test.rs[5,5],
                  heavy.static.test.rs[5,4]))



