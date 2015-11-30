install.packages(с("dplyr", "xts", "rugarch", "forecast", "car"),
                 repos = 'http://cran.us.r-project.org')

require(dplyr)
require(xts)
require(rugarch)
require(forecast)
require(car)


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
# one-step-ahead without re-estimation


heavy.static.test.sp <-
    data.frame(dtaTestSP[ ,"rcto"], dtaTestSP[ ,"rk"]) %>%
    setNames(c("rcto", "rk")) %>%
    mutate(garch = VarStaticForc(init.var = dtaEViewsSP$rcto %>% var(),
                                 const    = heavy.static.coef.sp[1,  "Coef"],
                                 coef     = heavy.static.coef.sp[3:4,"Coef"],
                                 coef.var = heavy.static.coef.sp[2,  "Coef"],
                                 data     = cbind(rk, rcto) %>% head(-1)),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch) %>% tail(-1)

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
           qlike = log(garch) - rk / garch) %>% tail(-1)


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
           qlike = log(garch) - rk / garch) %>% tail(-1)


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
           qlike = log(garch) - rk / garch) %>% tail(-1)


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
rgarch.static.forc.sp <-
    ugarchforecast(rgarch.specf.sp,
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
rgarch.static.forc.rs <-
    ugarchforecast(rgarch.specf.rs,
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

#---S&P 500---
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

garch.static.forc.sp <-
    ugarchforecast(garch.specf.sp,
                   n.ahead = 1,
                   n.roll = nrow(dtaTestSP[ ,"rcto"]) - 1,
                   data = dtaTestSP[ ,"rcto"],
                   out.sample = nrow(dtaTestSP[ ,"rcto"]) - 1) %>%
    sigma() %>% t() %>% as.data.frame() %>% head(-1) %>% setNames("sigma") %>%
    mutate(garch = sigma^2,
           rk    = tail(dtaTestSP[ ,"rk"], - 1),
           mse   = (rk - garch)^2,
           qlike = log(garch) - rk / garch)


#---RUSSEL---
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

garch.static.forc.rs <-
    ugarchforecast(garch.specf.rs,
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
## Some cleaning up
## ----------------------------------------------------------------------


clean.num.sp <- c(148, 149, 188:198, 247, 380)
clean.num.rs <- c(188:198, 225, 226, 227, 229, 230, 321, 322)

heavy.static.clean.sp  <- heavy.static.test.sp[-clean.num.sp, ]
mem.static.clean.sp    <- mem.static.test.sp[-clean.num.sp, ]
rgarch.static.clean.sp <- rgarch.static.forc.sp[-clean.num.sp, ]
garch.static.clean.sp  <- garch.static.forc.sp[-clean.num.sp, ]

heavy.static.clean.rs  <- heavy.static.test.rs[-clean.num.rs, ]
mem.static.clean.rs    <- mem.static.test.rs[-clean.num.rs, ]
rgarch.static.clean.rs <- rgarch.static.forc.rs[-clean.num.rs, ]
garch.static.clean.rs  <- garch.static.forc.rs[-clean.num.rs, ]


##------------------------------------------------------------------------
## Ststic one-step-ahead final comparison
## ----------------------------------------------------------------------

# Tips fo DM test:
# if DM<0 then loss(e1)<loss(e2), where e2 - bench
# if DM>0 then loss(e1)>loss(e2), where e2 - bench

SOSAFIC <- function(realized.measure, model.forecasts,
                        bench.forecast, model.names){

    # "LR" - linear regression coefficients,
    # "MZ" - Minzer-Zarnowiz test,
    # "LF" - MSE, QLIKE loss functions,
    # "DM" - Debold-Mariano test


    n.row = nrow(model.forecasts)
    n.col = ncol(model.forecasts)

    mse.model = matrix(ncol = n.col+1, nrow = n.row)
    qlike.model = matrix(ncol = n.col+1, nrow = n.row)
    loss = matrix(nrow = 2, ncol = n.col+1)

    
    colnames(loss) <- c(model.names, "bench")
    rownames(loss) <- c("MSE", "QLIKE")

    tmp.forecasts = cbind(model.forecasts, bench.forecast)

    ## LF for models
    for(i in 1:(n.col+1)){
        mse.model[,i] = (realized.measure - tmp.forecasts[,i])^2
        qlike.model[,i] = log(tmp.forecasts[,i]) +
                              realized.measure/tmp.forecasts[,i]

        loss["MSE", i] = mean(mse.model[,i])
        loss["QLIKE", i] = mean(qlike.model[,i])
    }

    ## DM-test
    dm = matrix(nrow = 2, ncol = n.col)
    colnames(dm) <- model.names
    rownames(dm) <- c("MSE", "QLIKE")
    for(i in 1:n.col){
        dm["MSE",i] = dm.test(e1 = mse.model[,i],
                              e2 = mse.model[,n.col+1],
                              #e2 = mse.bench,
                              alternative = "two.sided",
                              h = 1,
                              power = 1)["statistic"]%>%as.numeric()
        dm["QLIKE",i] = dm.test(e1 = qlike.model[,i],
                                e2 = qlike.model[,n.col+1],
                                #e2 = qlike.bench,
                                alternative = "two.sided",
                                h = 1,
                                power = 1)["statistic"]%>%as.numeric()
    }

    ## Mincer-Zarnowitz Linear Regression
    lr = matrix(nrow = 2, ncol = n.col+1)
    mz = matrix(nrow = 2, ncol = n.col+1)
    colnames(lr) <- c(model.names, "bench")
    rownames(lr) <- c("a0", "a1")
    colnames(mz) <- c(model.names, "bench")
    rownames(mz) <- c("ub", "mz") #ub - unbiasness test
    for(i in 1:(n.col+1)){
        lm.tmp = lm(realized.measure ~ tmp.forecasts[,i])
        names.tmp = names(lm.tmp$coefficients)

        lr[1, i] = lm.tmp$coefficients[1]
        lr[2, i] = lm.tmp$coefficients[2]

        mz[1, i] = linearHypothesis(
                       model = lm.tmp,
                       c(paste(names.tmp[1],'= 0')))$F[2]
        mz[2, i] = linearHypothesis(
                       model = lm.tmp,
                       c(paste(names.tmp[1],'= 0'),
                         paste(names.tmp[2],'= 1')))$F[2]
    }

    return(list("LR" = lr, "MZ" = mz, "LF" = loss, "DM" = dm))
}



comparison.sp <-
    SOSAFIC(realized.measure = heavy.static.clean.sp[,"rk"],
            model.forecasts = matrix(c(
                              heavy.static.clean.sp[,"garch"],
                              mem.static.clean.sp[,"garch"],
                              rgarch.static.clean.sp[,"garch"]), ncol=3),
            bench.forecast = garch.static.clean.sp[,"garch"],
            model.names = c("heavy", "mem", "rgarch"))


comparison.rs <-
    SOSAFIC(realized.measure = heavy.static.clean.rs[,"rk"],
            model.forecasts = matrix(c(
                              heavy.static.clean.rs[,"garch"],
                              mem.static.clean.rs[,"garch"],
                              rgarch.static.clean.rs[,"garch"]), ncol=3),
            bench.forecast = garch.static.clean.rs[,"garch"],
            model.names = c("heavy", "mem", "rgarch"))
