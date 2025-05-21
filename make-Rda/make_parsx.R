# Creates parameter objects

# Some temporary pars for testing
grp_parsx <- list(grps = c('m0', 'm1', 'm2','sr1'),
                 yield = c(default = 0.05, sr1 = 0.065),
                 xa_fresh = c(m0 = 0.0628, m1 = 0.0628, m2 = 0.0628, sr1 = 0.0628),
                 xa_init = c(all = 0.0628),
                 decay_rate = c(all = 0.02),
                 ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                 qhat_opt = c(m0 = 0.693, m1 = 0.407, m2 = 1.65, m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8.95),
                 T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                 T_min = c(m0 = 0, m1 = 6.41, m2 = 6.41, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                 T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
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
                  decay_rate_xa = 0.02)

man_parsx <- list(conc_fresh = list(sulfide = 0.01, 
				    sulfate = 0.2, 
				    TAN = 0.0, 
                                    VFA = 1.7, 
				    xa_aer = 0, 
				    xa_bac = 0, 
				    xa_dead = 0, 
				    VSd = 0, 
				    ash = 15), 
                  pH = 7, dens = 1000)

chem_parsx <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561,
                                VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, 
				CO2_aer = 1/0.436, CO2_sr = 1/1.2, 
                                C_xa = 1/0.3753125))



save(grp_parsx, file = '../data/grp_parsx.rda')
save(mic_parsx, file = '../data/mic_parsx.rda')
save(man_parsx, file = '../data/man_parsx.rda')
save(chem_parsx, file = '../data/man_parsx.rda')

