install.packages("xts")
library(xts)

# Data from 'Oxford-Man Institute of Quantitative Finance'. This data based on underlying
# high frequency data, which obtained through 'Reuters DataScope Tick History'.
# 
# # Rules of cleaning up of this data: 
# # http://realized.oxford-man.ox.ac.uk/documentation/data-cleaning

if(!dir.exists("../data")){dir.create("../data")}
fURL <- "http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip"
download.file(fURL, "../data/rowOxfordData.zip");
unzip("../data/rowOxfordData.zip", exdir = "../data")

rawOxfordData <- xts(read.zoo("../data/OxfordManRealizedVolatilityIndices.csv",  
                              format = "%Y%m%d", index.column = 1, 
                              header = T, sep= ",", skip = 2))
                             

