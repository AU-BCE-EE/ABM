# Check COD balance

checkCOD <- function(dat, 
                     grps,
                     subs,
                     COD_conv,
                     stoich,
                     rtol = 0.001
                    ) {

  first <- unlist(dat[1, ])
  last <- unlist(dat[nrow(dat), ])

  CODin <- sum(last[['COD_load']], first[grps], first[subs] * stoich['VFA', subs], first[['VFA']])
  CODeff <- sum(last[paste0(subs, '_eff')] * stoich['VFA', subs]) + last[['VFA_eff']]
  CODemis <- last[['CH4_emis_cum']] * COD_conv[['CH4']]
  CODrem <- sum(last[grps], last[subs] * stoich['VFA', subs], last[['VFA']])
  bal <- CODin - CODeff - CODemis - CODrem
  rbal <- bal / CODin

  if (abs(rbal) > rtol) {
    warning('COD balance is off by ', signif(100 * rbal, 2), '%')
    return(invisible(rbal))
  } 

  return(invisible(rbal))

}
