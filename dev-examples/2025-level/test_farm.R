rm(list = ls())

library(readxl)
library(dplyr)
library(ABM)

#dat <- '../inputs/saman_paper/Sow lactating average diet.xlsm'
#dat <- '../inputs/saman_paper/Sow gestation average diet.xlsm'
#dat <- '../inputs/saman_paper/Growing pig average diet.xlsm'
#dat <- '../inputs/saman_paper/Growing pig 5% sugar beet.xlsm'
#dat <- '../inputs/vanderzaag_2019/vanderzaag_2019.xlsm'
dat <- '../inputs/main_sheet.xlsm'
#dat <- 'C:/Users/au277187/OneDrive - Aarhus universitet/Documents/GitHub/ABM_parm/inputs/kari_2018.xlsm'

ff <- list.files('../R', full.names = T)
for(i in ff) source(i)

#ff <- list.files('../../ABM/R', full.names = T)
#for(i in ff) source(i)


#
#debug(abm_farm)
#debug(get_farm)
#debug(rates)
#debug(slurry_app_mod)
#debug(calcNorm)
out <- abm_farm(dat = dat, storage_mode = TRUE, years = 2)



