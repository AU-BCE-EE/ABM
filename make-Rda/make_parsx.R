# Creates parameter objects

# Some temporary pars for testing
grp_parsx <- list(grps = c('m0', 'm1', 'm2','sr1'),
                 yield = c(default = 0.05, sr1 = 0.065),
                 xa_fresh = c(m0 = 0.0628, m1 = 0.0628, m2 = 0.0628, sr1 = 0.0628),
                 xa_init = c(all = 0.0628),
                 dd_rate = c(all = 0.02),
                 ks = c(default = 1.153337, sr1 = 0.461335),
                 qhat_opt = c(m0 = 1, m1 = 1, m2 = 2, sr1 = 8.95),
                 T_opt = c(m0 = 18, m1 = 18, m2 = 28, sr1 = 43.75),
                 T_min = c(m0 = 0, m1 = 6.41, m2 = 6.41, sr1 = 0),
                 T_max = c(m0 = 25, m1 = 25, m2 = 38, sr1 = 51.25),
                 ki_NH3_min = c(all = 0.015),
                 ki_NH3_max = c(all = 0.3),
                 ki_NH4_min = c(all = 2.7),
                 ki_NH4_max = c(all = 7),
                 ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                 ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                 ki_H2S_min = c(default = 0.08),
                 IC50_low = c(default = 0.20854, sr1 = 0.2772),
                 ki_HAC = c(default = 0.31),
                 pH_UL = c(default = 8.0),
                 pH_LL = c(default = 6.5, sr1 = 5.5))

mic_parsx <- list(ks_SO4 = 0.00694, 
                  km_urea = 0.913,
                  dd_rate_xa = 0.02,
                  T_opt_hyd = 60,
                  T_min_hyd = 0,
                  T_max_hyd = 90,
                  hydrol_opt = 0.1)

man_parsx <- list(conc_fresh = c(sulfide = 0.01, 
				                         sulfate = 0.2, 
				                         TAN = 0.0, 
                                 VFA = 2, 
				                         VSd = 50, 
				                         ash = 15), 
                  pH = 7, dens = 1000)

chem_parsx <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561,
                                VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, 
                                CO2_aer = 1/0.436, CO2_sr = 1/1.2, 
                                C_xa = 1/0.3753125))



save(grp_parsx, file = '../data/grp_parsx.rda')
save(mic_parsx, file = '../data/mic_parsx.rda')
save(man_parsx, file = '../data/man_parsx.rda')
save(chem_parsx, file = '../data/chem_parsx.rda')

