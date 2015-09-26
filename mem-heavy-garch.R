install.packages(с("dplyr", "xts", "rugarch"))
require(dplyr)
require(xts)
require(rugarch)

##---------------------------------------------------------------------------
## Import HEAVY and MEM estimators from EViews through CSV
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
    data.frame(dtaTestSP[ ,"rcto"], dtaTestSP[ ,"rk"]) %>%
    setNames(c("rcto", "rk")) %>%
    mutate(garch = VarStaticForc(init.var = dtaEViewsSP$rcto %>% var(),
                                 const    = heavy.static.coef.sp[1,  "Coef"],
                                 coef     = heavy.static.coef.sp[3:4,"Coef"],
                                 coef.var = heavy.static.coef.sp[2,  "Coef"],
                                 data     = cbind(rk, rcto) %>% head(-1)),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


heavy.static.test.rs <-
    data.frame(dtaTestRS[ ,"rcto"], dtaTestRS[ ,"rk"]) %>%
    setNames(c("rcto", "rk")) %>%
    mutate(ind.rk = ifelse(rcto < 0, 1, 0) * rk,
           garch  = VarStaticForc(init.var = dtaEViewsRUS$rcto %>% var(),
                                  const    = heavy.static.coef.rs[1,  "Coef"],
                                  coef     = heavy.static.coef.rs[3:4,"Coef"],
                                  coef.var = heavy.static.coef.rs[2,  "Coef"],
                                  data     = cbind(ind.rk, rcto) %>% head(-1)),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


##----------------------------------------------------------------------------
## MEM forecast
##----------------------------------------------------------------------------
# 1-step-ahead without re-estimation


mem.static.test.sp <-
    data.frame(dtaTestSP[ ,"rcto"],
               dtaTestSP[ ,"rk"],
               dtaTestSP[ ,"hilo"]) %>%
    setNames(c("rcto", "rk", "hilo")) %>%
    mutate(hilo2 = hilo^2,
           garch = VarStaticForc(init.var = dtaEViewsSP$rcto %>% var(),
                                 const    = mem.static.coef.sp[1,  "Coef"],
                                 coef     = mem.static.coef.sp[3:5,"Coef"],
                                 coef.var = mem.static.coef.sp[2,  "Coef"],
                                 data     = cbind(rk, hilo2, rcto)%>% head(-1)),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


mem.static.test.rs <-
    data.frame(dtaTestRS[ ,"rcto"], dtaTestRS[ ,"rk"]) %>%
    setNames(c("rcto", "rk")) %>%
    mutate(ind = ifelse(rcto < 0, 1, 0),
           ind_rcto2 = ind * rcto^2,
           ind_rk = ind * rk,
           garch = VarStaticForc(init.var = dtaEViewsRUS$rcto %>% var(),
                                 const    = mem.static.coef.rs[1,  "Coef"],
                                 coef     = mem.static.coef.rs[3:5,"Coef"],
                                 coef.var = mem.static.coef.rs[2,  "Coef"],
                                 data     = cbind(ind_rcto2, rcto, ind_rk) %>%
                                            head(-1)),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)

##---------------------------------------------------------------------------
## Realized GARCH forecast
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
    ugarchfit(data = dtaTrainSP[ ,"rcto"],
              spec = rgarch.spec.sp,
              solver = "hybrid",
              realizedVol = dtaTrainSP[ ,"rk"] %>% sqrt())

rgarch.static.coef.rs <-
    ugarchfit(data = dtaTrainRS[ ,"rcto"],
              spec = rgarch.spec.rs,
              solver = "hybrid",
              realizedVol = dtaTrainRS[ ,"rk"] %>% sqrt())


#forecast
rgarch.specf.sp <- rgarch.spec.sp
setfixed(rgarch.specf.sp) <- as.list(coef(rgarch.static.coef.sp))
rgarch.static.forc.sp <- ugarchforecast(
    rgarch.specf.sp,
    n.ahead = 1,
    n.roll = nrow(dtaTestSP[ ,"rcto"]) - 1,
    data = dtaTestSP[ ,"rcto"],
    out.sample = nrow(dtaTestSP[ ,"rcto"]) - 1,
    realizedVol = dtaTestSP[ ,"rk"] %>% sqrt()) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk    = tail(dtaTestSP[ ,"rk"], -1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


rgarch.specf.rs <- rgarch.spec.rs
setfixed(rgarch.specf.rs) <- as.list(coef(rgarch.static.coef.rs))
rgarch.static.forc.rs <- ugarchforecast(
    rgarch.specf.rs,
    n.ahead = 1,
    n.roll = nrow(dtaTestRS[ ,"rcto"]) - 1,
    data = dtaTestRS[ ,"rcto"],
    out.sample = nrow(dtaTestRS[ ,"rcto"]) - 1,
    realizedVol = dtaTestRS[ ,"rk"] %>%  sqrt()) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk    = tail(dtaTestRS[ ,"rk"], -1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)

##---------------------------------------------------------------------------
## Simple GARCH forecast
##---------------------------------------------------------------------------
# 1-step-ahead without re-estimation

#---S&P 500
garch.spec.sp <-
    ugarchspec(variance.model = list(model = "sGARCH",
                                     garchOrder = c(1,1)),
               mean.model = list(armaOrder = c(0,0),
                                 include.mean = FALSE),
               distribution.model = "norm"
               #,fixed.pars = list("omega" = 0, "alpha1" = 0, "beta1" = 0)
              )

garch.static.coef.sp <-
    ugarchfit(data = dtaTrainSP[ ,"rcto"],
              spec = garch.spec.sp)

garch.specf.sp <- garch.spec.sp
setfixed(garch.specf.sp) <- as.list(coef(garch.static.coef.sp))

garch.static.forc.sp <- ugarchforecast(
    garch.specf.sp,
    n.ahead = 1,
    n.roll = nrow(dtaTestSP[ ,"rcto"]) - 1,
    data = dtaTestSP[ ,"rcto"],
    out.sample = nrow(dtaTestSP[ ,"rcto"]) - 1) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk    = tail(dtaTestSP[ ,"rk"], - 1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


#---RUSSEL
garch.spec.rs <-
    ugarchspec(variance.model = list(model = "sGARCH",
                                     garchOrder = c(1,1)),
               mean.model = list(armaOrder = c(0,0),
                                 include.mean = FALSE),
               distribution.model = "norm"
               ,fixed.pars = list("omega" = 0)
              )

garch.static.coef.rs <-
    ugarchfit(data = dtaTrainRS[ ,"rcto"],
              spec = garch.spec.rs)

garch.specf.rs <- garch.spec.rs
setfixed(garch.specf.rs) <- as.list(coef(garch.static.coef.rs))

garch.static.forc.rs <- ugarchforecast(
    garch.specf.rs,
    n.ahead = 1,
    n.roll = nrow(dtaTestRS[ ,"rcto"]) - 1,
    data = dtaTestRS[ ,"rcto"],
    out.sample = nrow(dtaTestRS[ ,"rcto"]) - 1,) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk    = tail(dtaTestRS[ ,"rk"], -1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


##------------------------------------------------------------------------
## Forecast function
## ----------------------------------------------------------------------

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

