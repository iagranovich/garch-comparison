install.packages(c("xts", "dplyr"))
library(xts)
library(dplyr)

# The dataset comes from 'Oxford-Man Institute of Quantitative Finance'
#
# # This data based on underlying high frequency data, which
# # obtained through 'Reuters DataScope Tick History'
#
# # The dataset contains log returns and realized measures
#
# # Rules of cleaning up of this data:
# # http://realized.oxford-man.ox.ac.uk/documentation/data-cleaning
#
# # Version of this dataset is 0.2

##------------------------------------------------------------------
## Dowloading and unzipping data
##------------------------------------------------------------------

if(!dir.exists("../data")){
    dir.create("../data")
}

fUrl <- "http://realized.oxford-man.ox.ac.uk/media/
          1366/oxfordmanrealizedvolatilityindices.zip"

download.file(fUrl, "../data/rowOxfordData.zip");

unzip("../data/rowOxfordData.zip", exdir = "../data")


##------------------------------------------------------------------
## Import row data from .CSV
##------------------------------------------------------------------

raw.file <- "../data/OxfordManRealizedVolatilityIndices.csv"


dtaOxfordRaw <-
    xts(read.zoo(file   = raw.file,
                 format = "%Y%m%d",
                 index.column = 1,
                 header = TRUE,
                 sep    = ",",
                 skip   = 2))


dtaNameColumns <-
    read.table(file  = raw.file,
               sep   = ",",
               nrows = 3) %>%
               t %>% tail(-1)


dtaShortName <- dtaNameColumns[grep("r$", dtaNameColumns[ ,3]), c(1,3)]


##--------------------------------------------------------------
## Extracting in-sample and out-of-sample data
##--------------------------------------------------------------

dtaTrain <- dtaOxfordRaw["2004-01-01/2007-12-31"]   #for estimation
dtaTest  <- dtaOxfordRaw["2008-01-01/2009-12-31"]   #for testing

##--------------------------------------------------------------
## Export data for calculating coeff by EViews
##--------------------------------------------------------------

# @param rcto  -- close-to-open returns
# @param rcto2 -- square close-to-open returns
# @param hilo  -- high-to-low returns
# @param hilo2 -- square high-to-low returns
# @param rk    -- realized kernel

dtaEViewsSP <-
    dtaTrain %>%
    as.data.frame() %>%
    select(open  = SPX2.openprice,
           close = SPX2.closeprice,
           rk    = SPX2.rk,
           hilo  = SPX2.highlow) %>%
    na.omit() %>%
    mutate(rcto  = log(close/open),
           rcto2 = rcto^2,
           hilo2 = hilo^2,
           ind   = ifelse(rcto < 0, 1, 0))


dtaEViewsRUS <-
    dtaTrain %>%
    as.data.frame() %>%
    select(open  = RUT2.openprice,
           close = RUT2.closeprice,
           rk    = RUT2.rk,
           hilo  = RUT2.highlow) %>%
    na.omit() %>%
    mutate(rcto  = log(close/open),
           rcto2 = rcto^2,
           hilo2 = hilo^2,
           ind   = ifelse(rcto < 0, 1, 0))


write.table(x = dtaEViewsSP,
            file = "../data/eviews_2008/sp2008_dtatrain.csv",
            sep  = ",",
            row.names = F)

write.table(x = dtaEViewsRUS,
            file = "../data/eviews_2008/rus2008_dtatrain.csv",
            sep  = ",",
            row.names = F)
