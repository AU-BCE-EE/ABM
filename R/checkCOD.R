# Check COD balance

checkCOD <- function(dat, 
                     grps = c('m0', 'm1', 'm2', 'sr1'),
                     subs = c('VSd'),
                     COD_conv = c(CH4 = 1/0.251, 
                                  xa = 1/0.707, 
                                  VFA = 1/0.938, 
                                  VSd = 1/0.69),
                     rtol = 0.001
                    ) {

  first <- unlist(dat[1, ])
  last <- unlist(dat[nrow(dat), ])

  CODin <- sum(last[['COD_load']], first[grps], first[subs], first[['VFA']])
  CODeff <- sum(last[paste0(subs, '_eff')]) + last[['VFA_eff']]
  CODemis <- last[['CH4_emis_cum']] * COD_conv[['CH4']]
  CODrem <- sum(last[grps], last[subs], last[['VFA']])
  bal <- CODin - CODeff - CODemis - CODrem
  rbal <- bal / CODin

  if (abs(rbal) > rtol) {
    warning('COD balance is off by ', signif(100 * rbal, 2), '%')
    return(invisible(rbal))
  } 

  return(invisible(rbal))

}
