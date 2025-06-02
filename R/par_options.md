
mng_pars = list(slurry_prod_rate = 10000, 
                slurry_mass = 1000,     
                storage_depth = 2,     
                resid_depth = 0.1,      
                area = 100,              
                empty_int = 100,          
                temp_C = 20,
                wash_water = 0,            
                wash_int = NA,
                rest_d = 0,
                resid_enrich = 1)

sub_pars <- list(subs = c('VSd'),
                 T_opt_hyd = c(VSd = 60),
                 T_min_hyd = c(VSd = 0),
                 T_max_hyd = c(VSd = 90),
                 hydrol_opt = c(VSd = 0.1),
                 sub_fresh = c(VSd = 50),
                 sub_init = c(VSd = 50))

grp_pars <- list(grps = c('m0', 'm1', 'm2','sr1'),
                 yield = c(default = 0.05, sr1 = 0.065),
                 xa_fresh = c(all = 0.05),
                 xa_init = c(all = 0.05),
                 dd_rate = c(all = 0.02),
                 ks = c(default = 1, sr1 = 0.5),
                 qhat_opt = c(m0 = 1, m1 = 1, m2 = 2, sr1 = 9),
                 T_opt = c(m0 = 18, m1 = 18, m2 = 28, sr1 = 44),
                 T_min = c(m0 = 0, m1 = 6.41, m2 = 6.41, sr1 = 0),
                 T_max = c(m0 = 25, m1 = 25, m2 = 38, sr1 = 51))

mic_pars <- list(ks_SO4 = 0.00694, 
                 km_urea = 0.913,
                 dd_rate_xa = 0.02)

man_pars <- list(conc_fresh = c(VFA = 2), 
                 pH = 7, dens = 1000)

chem_pars <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561,
                               VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, 
                               CO2_aer = 1/0.436, CO2_sr = 1/1.2, 
                               C_xa = 1/0.3753125))


Or should all components go here?

man_pars <- list(conc_fresh = c(VSd = ,
                                m1 = ,
                                ...
                                VFA = 2), 
                 pH = 7, dens = 1000)


