install.packages(с("dplyr", "xts", "rugarch"))
require(dplyr)
require(xts)
require(rugarch)

##---------------------------------------------------------------------------
## Import HEAVY estimators from EViews through CSV
##---------------------------------------------------------------------------

heavy.static.coef.sp <- read.csv("./coef_2008/heavy_coef_sp.csv")
heavy.static.coef.rs <- read.csv("./coef_2008/heavy_coef_rs.csv")
mem.static.coef.sp <- read.csv("./coef_2008/mem_coef_sp.csv")
mem.static.coef.rs <- read.csv("./coef_2008/mem_coef_rs.csv")

##---------------------------------------------------------------------------
## HAVY forecast
##---------------------------------------------------------------------------
# 1-step-ahead without re-estimation

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

##----------------------------------------------------------------------------
## MEM forecast
##----------------------------------------------------------------------------
# 1-step-ahead without re-estimation

mem.static.test.sp <-
    cbind(dtaTest[ ,"SPX2.openprice"],
          dtaTest[ ,"SPX2.closeprice"],
          dtaTest[ ,"SPX2.highlow"],
          dtaTest[ ,"SPX2.rk"]) %>%
    na.omit() %>%
    data.frame() %>%
    mutate(rcto  = log(SPX2.closeprice / SPX2.openprice),
           hilo2 = SPX2.highlow^2,
           garch = VarStaticForc(init.var = dtaEViewsSP$rcto %>% var(),
                                 const    = mem.static.coef.sp[1,  "Coef"],
                                 coef     = mem.static.coef.sp[3:5,"Coef"],
                                 coef.var = mem.static.coef.sp[2,  "Coef"],
                                 data     = cbind(SPX2.rk, hilo2, rcto) %>%
                                            head(-1)),
           mse   = (SPX2.rk - garch)^2,
           qlike = log(garch) - SPX2.rk / garch)

mem.static.test.rs <-
    cbind(dtaTest[ ,"RUT2.openprice"],
          dtaTest[ ,"RUT2.closeprice"],
          dtaTest[ ,"RUT2.rk"]) %>%
    na.omit() %>%
    data.frame() %>%
    mutate(rcto  = log(RUT2.closeprice / RUT2.openprice),
           ind   = ifelse(rcto < 0, 1, 0),
           ind_rcto2 = ind * rcto^2,
           ind_rk = ind * RUT2.rk,
           garch = VarStaticForc(init.var = dtaEViewsRUS$rcto %>% var(),
                                 const    = mem.static.coef.rs[1,  "Coef"],
                                 coef     = mem.static.coef.rs[3:5,"Coef"],
                                 coef.var = mem.static.coef.rs[2,  "Coef"],
                                 data     = cbind(ind_rcto2, rcto, ind_rk) %>%
                                            head(-1)),
           mse   = (RUT2.rk - garch)^2,
           qlike = log(garch) - RUT2.rk / garch)

##---------------------------------------------------------------------------
## RGARCH forecast
##---------------------------------------------------------------------------
## 1-step-ahead without re-estimation

#spec
rgarch.spec.sp <-
    ugarchspec(variance.model = list(model = "realGARCH",
                                     garchOrder = c(1,1)),
               mean.model = list(armaOrder = c(0,0),
                                 include.mean = FALSE),
               distribution.model = "norm"
               ,fixed.pars = list("xi" = 0)
              )


rgarch.spec.rs <-
    ugarchspec(variance.model = list(model = "realGARCH",
                                     garchOrder = c(1,1)),
               mean.model = list(armaOrder = c(0,0),
                                 include.mean = FALSE),
               distribution.model = "norm"
               #,fixed.pars = list("omega" = 0)
              )
#fit
rgarch.static.coef.sp <-
    ugarchfit(data = log(dtaTrain[ ,"SPX2.closeprice"]/
                         dtaTrain[ ,"SPX2.openprice"] %>%
                         na.omit()),
              spec = rgarch.spec.sp,
              solver = "hybrid",
              realizedVol = dtaTrain[ ,"SPX2.rk"] %>%
                            na.omit() %>% sqrt())

rgarch.static.coef.rs <-
    ugarchfit(data = log(dtaTrain[ ,"RUT2.closeprice"]/
                         dtaTrain[ ,"RUT2.openprice"] %>%
                         na.omit()),
              spec = rgarch.spec.rs,
              solver = "hybrid",
              realizedVol = dtaTrain[ ,"RUT2.rk"] %>%
                            na.omit() %>% sqrt())

#forecast
rgarch.specf.sp <- rgarch.spec.sp
setfixed(rgarch.specf.sp) <- as.list(coef(rgarch.static.coef.sp))
rgarch.static.forc.sp <- ugarchforecast(
    rgarch.specf.sp,
    n.ahead = 1,
    n.roll = nrow(na.omit(dtaTest[ ,"SPX2.closeprice"]))-1,
    data = log(dtaTest[ ,"SPX2.closeprice"]/
               dtaTest[ ,"SPX2.openprice"] %>%
               na.omit()),
    out.sample = nrow(na.omit(dtaTest[ ,"SPX2.closeprice"]))-1,
    realizedVol = (dtaTest[ ,"SPX2.rk"] %>%
                   na.omit() %>% sqrt())) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk = tail(na.omit(dtaTest[ ,"SPX2.rk"]), -1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


rgarch.specf.rs <- rgarch.spec.rs
setfixed(rgarch.specf.rs) <- as.list(coef(rgarch.static.coef.rs))
rgarch.static.forc.rs <- ugarchforecast(
    rgarch.specf.rs,
    n.ahead = 1,
    n.roll = nrow(na.omit(dtaTest[ ,"RUT2.closeprice"]))-1,
    data = log(dtaTest[ ,"RUT2.closeprice"]/
               dtaTest[ ,"RUT2.openprice"] %>%
               na.omit()),
    out.sample = nrow(na.omit(dtaTest[ ,"RUT2.closeprice"]))-1,
    realizedVol = (dtaTest[ ,"RUT2.rk"] %>%
                   na.omit() %>% sqrt())) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk = tail(na.omit(dtaTest[ ,"RUT2.rk"]), -1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)

##---------------------------------------------------------------------------
## Simple GARCH forecast
##---------------------------------------------------------------------------
# 1-step-ahead without re-estimation




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


CheckStaticForc(mem.static.test.sp[494,7],
                mem.static.test.sp[493,7],
                mem.static.coef.sp[1,"Coef"],
                mem.static.coef.sp[3:5,"Coef"],
                mem.static.coef.sp[2,"Coef"],
                c(mem.static.test.sp[493,4],
                  mem.static.test.sp[493,6],
                  mem.static.test.sp[493,5]))

CheckStaticForc(mem.static.test.rs[7,8],
                mem.static.test.rs[6,8],
                mem.static.coef.rs[1,"Coef"],
                mem.static.coef.rs[3:5,"Coef"],
                mem.static.coef.rs[2,"Coef"],
                c(mem.static.test.rs[6,6],
                  mem.static.test.rs[6,4],
                  mem.static.test.rs[6,7]))

