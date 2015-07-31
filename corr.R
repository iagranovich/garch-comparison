# corr <- dataout %>% as.data.frame %>% 
#                    select(ends_with(".r")) %>% 
#                    cor(use = "pairwise.complete.obs")


## 
## Correlations with significance levels 
##
install.packages("Hmisc")
library(Hmisc)

corr <- dtaTest %>% as.data.frame %>% 
                    select(ends_with(".r")) %>%
                    as.matrix %>%
                    rcorr(type = "pearson")
