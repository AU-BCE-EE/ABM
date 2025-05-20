hard_pars <- function(pars){
  
  pars$g_NH4 <- 0.7
  pars$temp_standard <- 298
  pars$temp_zero <- 273
  pars$pH_floor <- 7
  pars$i_meth <- grep('^[mp]', names(pars$qhat_opt))
  pars$i_sr <- grep('^sr', names(pars$qhat_opt))
  pars$n_mic <- length(pars$qhat_opt)

  pars$pH_inhib  <- 0 * pars$pH_LL + 1 
  pars$NH3_inhib <- 0 * pars$pH_LL + 1
  pars$NH4_inhib <- 0 * pars$pH_LL + 1
  pars$HAC_inhib <- 0 * pars$pH_LL + 1
  pars$H2S_inhib <- 0 * pars$pH_LL + 1
  
  return(pars)
}
