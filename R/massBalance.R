# COD, C, and N mass balance
# Not exported, called by abm() and abm_*()

# Molar masses used to convert CH4 and CO2 emission (g) to C (g)
C_per_CH4 <- 12.011 / 16.043
C_per_CO2 <- 12.011 / 44.01

# COD, C, and N content (g) of slurry components
# x: named numeric vector or data frame with state variables (g), e.g., storage or effluent
# cc: COD_conv parameter vector
# mic_names: names of microbial biomass elements/columns in x
# N_VSd: N in VSd (g N/g COD)
# Microbial biomass and xa_dead do not contain N here, consistent with N_load_cum
elementContent <- function(x, cc, mic_names, N_VSd) {
  mic <- if (is.data.frame(x)) rowSums(x[, mic_names, drop = FALSE]) else sum(x[mic_names])
  COD <- mic + x[['xa_dead']] + x[['RFd']] + x[['iNDF']] + x[['VSd']] + x[['starch']] +
         x[['CPs']] + x[['CPf']] + x[['Cfat']] + x[['VFA']]
  C <- (mic + x[['xa_dead']]) / cc[['C_xa']] + x[['RFd']] / cc[['C_RFd']] + x[['iNDF']] / cc[['C_iNDF']] +
       x[['VSd']] / cc[['C_VSd']] + x[['starch']] / cc[['C_starch']] + (x[['CPs']] + x[['CPf']]) / cc[['C_CP']] +
       x[['Cfat']] / cc[['C_Cfat']] + x[['VFA_C']] + x[['urea']] / cc[['C_N_urea']]
  N <- (x[['CPs']] + x[['CPf']]) / cc[['CP_N']] + x[['VSd']] * N_VSd + x[['urea']] + x[['TAN']]
  return(list(COD = unname(COD), C = unname(C), N = unname(N)))
}

# N in VSd (g N/g COD), from formula C15.89H26.3O9.2N (1 mol N per mol, see mol_COD in hard_pars())
NVSd <- function(pars) 14.007 / pars$mol_COD[['OM']]

# Add COD, C, and N removed in effluent (y.eff, output from emptyStore()) to cumulative state variables in y
# Removal is tracked in state variables (not only in _eff columns) so the balance is correct in every output row,
# independent of which row effluent is reported in
addEffCum <- function(y, y.eff, pars) {
  e <- y.eff
  names(e) <- sub('_eff$', '', names(e))
  # Microbial biomass is always first in the state variable vector
  el <- elementContent(e, pars$COD_conv, names(e)[seq_len(pars$n_mic)], NVSd(pars))
  y['COD_eff_cum'] <- y['COD_eff_cum'] + el$COD
  y['C_eff_cum'] <- y['C_eff_cum'] + el$C
  y['N_eff_cum'] <- y['N_eff_cum'] + el$N
  return(y)
}

# Mass balance for COD, C, and N for each row of output dat (after renaming . to _)
# For each element:
#   in       = initial storage + cumulative load
#   stored   = current storage
#   removed  = cumulative removal in effluent
#   emitted  = cumulative conversion/emission (COD: CH4, respiration, sulfate reduction; C: CH4, CO2; N: NH3, N2O)
#   gap      = cumulative net loss due to process stoichiometry (sum of *_gap_*_cum, positive = destroyed, negative = created)
#   resid    = in - stored - removed - emitted - gap (bookkeeping error, should be ~0)
#   rel_resid = resid / in
# All cumulative variables are taken relative to the first row, so a run can start from any state
massBalance <- function(dat, pars) {
  cc <- pars$COD_conv
  el <- elementContent(dat, cc, pars$grps, NVSd(pars))
  d0 <- function(x) x - x[1]
  gapSum <- function(el) {
    g <- dat[, grepl(paste0('^', el, '_gap_.+_cum$'), names(dat)), drop = FALSE]
    d0(rowSums(g))
  }

  bal <- list(
    COD = list(stored = el$COD,
               load = d0(dat$COD_load_cum),
               removed = d0(dat$COD_eff_cum),
               emitted = d0(dat$COD_conv_cum_meth) + d0(dat$COD_conv_cum_respir) + d0(dat$COD_conv_cum_sr),
               gap = gapSum('COD')),
    C = list(stored = el$C,
             load = d0(dat$C_load_cum),
             removed = d0(dat$C_eff_cum),
             emitted = d0(dat$CH4_emis_cum) * C_per_CH4 + d0(dat$CO2_emis_cum) * C_per_CO2,
             gap = gapSum('C')),
    N = list(stored = el$N,
             load = d0(dat$N_load_cum),
             removed = d0(dat$N_eff_cum),
             emitted = d0(dat$NH3_emis_cum) + d0(dat$N2O_emis_cum),
             gap = gapSum('N'))
  )

  out <- list()
  for (e in names(bal)) {
    b <- bal[[e]]
    b_in <- b$stored[1] + b$load
    resid <- b_in - b$stored - b$removed - b$emitted - b$gap
    out[[paste0(e, '_bal_in')]] <- b_in
    out[[paste0(e, '_bal_stored')]] <- b$stored
    out[[paste0(e, '_bal_removed')]] <- b$removed
    out[[paste0(e, '_bal_emitted')]] <- b$emitted
    out[[paste0(e, '_bal_gap')]] <- b$gap
    out[[paste0(e, '_bal_resid')]] <- resid
    out[[paste0(e, '_bal_rel_resid')]] <- resid / pmax(b_in, .Machine$double.eps)
  }

  return(as.data.frame(out))
}
