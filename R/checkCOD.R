# Check COD balance

dat <- out0
head(dat)
dat
checkCOD <- function(dat, 
                     grps = c('m0', 'm1', 'm2', 'sr1'),
                     subs = c('VSd'),
                     COD_conv = c(CH4 = 1/0.251, 
                                  xa = 1/0.707, 
                                  VFA = 1/0.938, 
                                  VSd = 1/0.69),
                     rtol = 0.001
                    ) {

  first <- dat[1, ]
  last <- dat[nrow(dat), ]

  CODin <- last$COD_load + first[[subs]] + first[['VFA']]
  #CODeff <- last[[paste0(subs, '_eff')]] + last[['VFA']]
  CODemis <- last[['CH4_emis_cum']] / COD_conv[['CH4']]
  CODrem <- sum(last[[subs]], last[['VFA']])
  bal <- CODin - CODemis - CODrem
  rbal <- bal / CODin
  #bal <- CODin - CODeff - CODemis - CODrem


}
