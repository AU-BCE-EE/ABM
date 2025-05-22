rates <- function(t, 
		              y, 
		              parms, 
		              temp_C_fun = temp_C_fun, 
		              pH_fun = pH_fun) {

    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Time-dependent values
  pH <- pH_fun(t + p$t_run)
  temp_C <- temp_C_fun(t + p$t_run)
  temp_K <- temp_C + 273.15

  # Hydrolysis rate . . . 
  alpha <- 0.1    

  # Microbial substrate utilization rate (vectorized calculation)
  qhat <- CTM_cpp(temp_K, p$T_opt, p$T_min, p$T_max, p$qhat_opt)
  names(qhat) <- names(p$qhat_opt)

  # VFA consumption rates (g/d)
  rut <- NA * qhat
  rut[p$i_meth] <- (qhat[p$i_meth] * y['VFA'] * y[p$i_meth]) / 
                   (p$ks[p$i_meth] * y['slurry_mass'] + y['VFA'])
  rut[p$i_sr] <- qhat[p$i_sr] * 0
 
  # Derivatives, all in g/d as COD except slurry_mass (kg/d as fresh slurry mass)
  # First element is xa, to avoid `xa.` prefix from names, simply omit xa = here
  ders <- c(p$yield * rut - p$dd_rate * y[1:p$n_mic],
            slurry_mass = p$slurry_prod_rate,
            VSd = p$slurry_prod_rate * p$conc_fresh[['VSd']] + 
                  sum(y[c(p$i_meth, p$i_sr)] * p$dd_rate) -
                  alpha * y['VSd'],
            VFA = p$slurry_prod_rate * p$conc_fresh[['VFA']] + 
                  alpha * y['VSd'] - 
                  sum(rut),
            CH4_emis_cum = sum(rut[p$i_meth]) / p$COD_conv[['CH4']],
            CO2_emis_cum = 0,
            slurry_load = p$slurry_prod_rate,
            COD_load = (p$conc_fresh[['VSd']] + p$conc_fresh[['VFA']]) * p$slurry_prod_rate)

  return(list(ders, 
         c(CH4_emis_rate = sum(rut[p$i_meth]) / p$COD_conv[['CH4']]), 
           temp_C = temp_C,
           pH = pH
          )
  )

}
