# Creates parameter objects

# best fit grp_pars2.0 based on multiple datasets using CP, RFd, starch and Cfat as inputs
grp_pars2.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                 yield = c(default = 0.05, sr1 = 0.065),
                 xa_fresh = c(default = 0.0628),
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
                 ki_NH4_max = c(all = 19.14),
                 ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                 ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                 ki_H2S_min = c(default = 0.08),
                 IC50_low = c(default = 0.20854, sr1 = 0.2772),
                 ki_HAC = c(default = 0.31),
                 pH_UL = c(default = 8.0),
                 pH_LL = c(default = 6.5, sr1 = 5.5))

mic_pars2.0 <- list(ks_SO4 = 0.00694, 
                    km_urea = 0.913,
                    decay_rate_xa = 0.02)

man_pars2.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.2, TAN = 0.0, starch = 5.25, 
                                  VFA = 1.7, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 27.6, CP = 21.1, RFd = 25.4, iNDF = 11.3, VSd = 0, 
                                  VSd_A = 55, VSnd_A = 23.5, ash = 15), pH = 7, dens = 1000)

wthr_pars2.0 <- list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 10)

chem_pars2.0 <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561, RFd = 1/0.8444792, iNDF = 1/0.8444792, starch = 1/0.8444792, 
                              Cfat = 1/0.3117844, CP = 1/0.6541602, VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, CO2_aer = 1/0.436, CO2_sr = 1/1.2, CO2_ureo = 1/1.57,
                              N_CP = 1/0.1014, C_xa = 1/0.3753125, C_RFd = 1/0.376, C_iNDF = 1/0.358, N_xa = 1/0.08754375,
                              C_starch = 1/0.377, C_Cfat = 1/0.265, C_CP = 1/0.359 , C_VFA = 1/0.374, C_VSd = 1/0.344, C_N_urea = 1/0.429,
                              frac_CP_xa = 0.835
                              ))


arrh_pars2.0 <- list(lnA = c(VSd_A = 31.3),
                     E_CH4 = c(VSd_A = 81000), 
                     A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                     E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                     R = 8.314,  
                     VS_CH4 = 6.67,
                     scale_alpha_opt = list(VSd = 1, notVSd = 1),
                     kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

save(grp_pars2.0, file = '../data/grp_pars2.0.rda')
save(mic_pars2.0, file = '../data/mic_pars2.0.rda')
save(man_pars2.0, file = '../data/man_pars2.0.rda')
save(wthr_pars2.0, file = '../data/wthr_pars2.0.rda')
save(chem_pars2.0, file = '../data/chem_pars2.0.rda')
save(arrh_pars2.0, file = '../data/arrh_pars2.0.rda')

# best fit grp_pars for Dalby et al. 2023, EST using VSd as input
grp_pars1.0 <- list(grps = c('m0','m1','m2', 'sr1'),
                               yield = c(default = 0.05, sr1 = 0.065),
                               xa_fresh = c(default = 0.0628),
                               xa_init = c(all = 0.0628),
                               decay_rate = c(all = 0.02),
                               ks_coefficient = c(default = 1.17071975112963, sr1 = 0.4682879),
                               qhat_opt = c(m0 = 0.4742862, m1 = 1.138287, m2 = 1.770668 , m3 = 7.2, m4 = 8, m5 = 8, sr1 = 2.529526),
                               T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                               T_min = c(m0 = 0, m1 = 10, m2 = 10, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                               T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                               ki_NH3_min = c(all = 0.015),
                               ki_NH3_max = c(all = 0.13),
                               ki_NH4_min = c(all = 2.7),
                               ki_NH4_max = c(all = 4.8),
                               ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                               ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                               ki_H2S_min = c(default = 0.08),
                               IC50_low = c(default = 0.20854, sr1 = 0.2772),
                               ki_HAC = c(default = 0.31),
                               pH_UL = c(default = 8.0),
                               pH_LL = c(default = 6.5, sr1 = 5.5))


man_pars1.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.2, TAN = 0.0, starch = 0, 
                                      VFA = 2, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 0, CP = 0, RFd = 0, iNDF = 0, VSd = 75, 
                                      VSd_A = 55, VSnd_A = 22, ash = 15), pH = 7, dens = 1000)

save(grp_pars1.0, file = '../data/grp_pars1.0.rda')
save(man_pars1.0, file = '../data/man_pars1.0.rda')

## set in progress 3.0
grp_pars_VS_pig3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                    yield = c(default = 0.05, sr1 = 0.065),
                    xa_fresh = c(default = 0.06),
                    xa_init = c(all = 0.06),
                    decay_rate = c(all = 0.02),
                    ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                    qhat_opt = c(m0 = 0.71216, m1 = 0.248, m2 = 1.7222, m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8.95),
                    T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                    T_min = c(m0 = 0, m1 = 12.59, m2 = 12.59, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                    T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                    ki_NH3_min = c(all = 0.015),
                    ki_NH3_max = c(all = 0.3576),
                    ki_NH4_min = c(all = 2.7),
                    ki_NH4_max = c(all = 19.16),
                    ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                    ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                    ki_H2S_min = c(default = 0.08),
                    IC50_low = c(default = 0.20854, sr1 = 0.2772),
                    ki_HAC = c(default = 0.31),
                    pH_UL = c(default = 8.0),
                    pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_pig3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                       yield = c(default = 0.05, sr1 = 0.065),
                       xa_fresh = c(default = 0.06),
                       xa_init = c(all = 0.06),
                       decay_rate = c(all = 0.02),
                       ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                       qhat_opt = c(m0 = 0.7152, m1 = 0.1797, m2 = 1.6975, m3 = 8, m4 = 8, m5 = 8, sr1 = 8.95),
                       T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                       T_min = c(m0 = 0, m1 = 12.54, m2 = 12.54, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                       T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                       ki_NH3_min = c(all = 0.015),
                       ki_NH3_max = c(all = 0.36187),
                       ki_NH4_min = c(all = 2.7),
                       ki_NH4_max = c(all = 20.05),
                       ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                       ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                       ki_H2S_min = c(default = 0.08),
                       IC50_low = c(default = 0.20854, sr1 = 0.2772),
                       ki_HAC = c(default = 0.31),
                       pH_UL = c(default = 8.0),
                       pH_LL = c(default = 6.5, sr1 = 5.5))




grp_pars_VS_cattle3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                              yield = c(default = 0.05, sr1 = 0.065),
                              xa_fresh = c(default = 0.3),
                              xa_init = c(all = 0.3),
                              decay_rate = c(all = 0.02),
                              ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                              qhat_opt = c(m0 = 0.815577, m1 = 0.01, m2 = 1.86466, m3 = 8.0, m4 = 8, m5 = 8, sr1 = 8.95),
                              T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                              T_min = c(m0 = 0, m1 = 5.04, m2 = 5.04, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                              T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                              ki_NH3_min = c(all = 0.015),
                              ki_NH3_max = c(all = 0.297),
                              ki_NH4_min = c(all = 2.7),
                              ki_NH4_max = c(all = 19.14),
                              ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                              ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                              ki_H2S_min = c(default = 0.08),
                              IC50_low = c(default = 0.20854, sr1 = 0.2772),
                              ki_HAC = c(default = 0.31),
                              pH_UL = c(default = 8.0),
                              pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_cattle3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                        yield = c(default = 0.05, sr1 = 0.065),
                        xa_fresh = c(default = 0.3),
                        xa_init = c(all = 0.3),
                        decay_rate = c(all = 0.02),
                        ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                        qhat_opt = c(m0 = 0.8607, m1 = 0.01, m2 = 1.513128, m3 = 8.0, m4 = 8, m5 = 8, sr1 = 8.95),
                        T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                        T_min = c(m0 = 0, m1 = 5.15, m2 = 5.15, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                        T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                        ki_NH3_min = c(all = 0.015),
                        ki_NH3_max = c(all = 0.3),
                        ki_NH4_min = c(all = 2.7),
                        ki_NH4_max = c(all = 19.14),
                        ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                        ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                        ki_H2S_min = c(default = 0.08),
                        IC50_low = c(default = 0.20854, sr1 = 0.2772),
                        ki_HAC = c(default = 0.31),
                        pH_UL = c(default = 8.0),
                        pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_VS_digestate3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                              yield = c(default = 0.05, sr1 = 0.065),
                              xa_fresh = c(default = 0.3),
                              xa_init = c(all = 0.3),
                              decay_rate = c(all = 0.02),
                              ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                              qhat_opt = c(m0 = 0.123, m1 = 0.01, m2 = 1.684, m3 = 8.0, m4 = 8, m5 = 8, sr1 = 8.95),
                              T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                              T_min = c(m0 = 0, m1 = 6.65, m2 = 6.65, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                              T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                              ki_NH3_min = c(all = 0.015),
                              ki_NH3_max = c(all = 0.3),
                              ki_NH4_min = c(all = 2.7),
                              ki_NH4_max = c(all = 19.14),
                              ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                              ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                              ki_H2S_min = c(default = 0.08),
                              IC50_low = c(default = 0.20854, sr1 = 0.2772),
                              ki_HAC = c(default = 0.31),
                              pH_UL = c(default = 8.0),
                              pH_LL = c(default = 6.5, sr1 = 5.5))


grp_pars_digestate3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                           yield = c(default = 0.05, sr1 = 0.065),
                           xa_fresh = c(default = 0.3),
                           xa_init = c(all = 0.3),
                           decay_rate = c(all = 0.02),
                           ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                           qhat_opt = c(m0 = 0.26, m1 = 0.01, m2 = 2.01, m3 = 8.0, m4 = 8, m5 = 8, sr1 = 8.95),
                           T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                           T_min = c(m0 = 0, m1 = 9.813, m2 = 9.81, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                           T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                           ki_NH3_min = c(all = 0.015),
                           ki_NH3_max = c(all = 0.3),
                           ki_NH4_min = c(all = 2.7),
                           ki_NH4_max = c(all = 19.14),
                           ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                           ki_H2S_int = c(default = 1.02, sr1 = 1.42),
                           ki_H2S_min = c(default = 0.08),
                           IC50_low = c(default = 0.20854, sr1 = 0.2772),
                           ki_HAC = c(default = 0.31),
                           pH_UL = c(default = 8.0),
                           pH_LL = c(default = 6.5, sr1 = 5.5))


mic_pars3.0 <- list(ks_SO4 = 0.00694,
                    km_urea = 0.913,
                    decay_rate_xa = 0.02
                    )
                    

man_pars_pig3.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.01, TAN = 0.0, starch = 5.25, 
                                      VFA = 1.7, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 27.6, CP = 21.1, RFd = 25.4, iNDF = 11.3, VSd = 0, 
                                      VSd_A = 55, VSnd_A = 23.5, ash = 15), pH = 7, dens = 1000)

man_pars_cattle3.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 2.116, sulfate = 0.01, TAN = 0.0, starch = 1.4864, 
                                          VFA = 2.54, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 16.04, CP = 30, RFd = 44.4, iNDF = 21.465, VSd = 0, 
                                          VSd_A = 32.86, VSnd_A = 45.4, ash = 15), pH = 7, dens = 1000)

man_pars_digestate3.0 <- list(conc_fresh = list(sulfide = 0, urea = 0, sulfate = 0, TAN = 3, starch = 0.29, 
                                             VFA = 0.5, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 3.18, CP = 5.95, RFd = 8.8, iNDF = 21.47, VSd = 0, 
                                             VSd_A = 12.6, VSnd_A = 18.13, ash = 15), pH = 8, dens = 1000)



man_pars_VS_pig3.0 <- man_pars_pig3.0
man_pars_VS_cattle3.0 <- man_pars_cattle3.0
man_pars_VS_digestate3.0 <- man_pars_digestate3.0

# pig
man_pars_VS_pig3.0$conc_fresh['VSd'] <- (man_pars_VS_pig3.0$conc_fresh[['starch']] + man_pars_VS_pig3.0$conc_fresh[['xa_dead']] +
man_pars_VS_pig3.0$conc_fresh[['Cfat']] + man_pars_VS_pig3.0$conc_fresh[['RFd']] + man_pars_VS_pig3.0$conc_fresh[['CP']] + 
man_pars_VS_pig3.0$conc_fresh[['iNDF']]) * 0.7

for (i in c('starch', 'xa_dead','Cfat','RFd','CP','iNDF')){
  man_pars_VS_pig3.0$conc_fresh[i] <- 0
}

# cattle
man_pars_VS_cattle3.0$conc_fresh['VSd'] <- (man_pars_VS_cattle3.0$conc_fresh[['starch']] + man_pars_VS_cattle3.0$conc_fresh[['xa_dead']] +
man_pars_VS_cattle3.0$conc_fresh[['Cfat']] + man_pars_VS_cattle3.0$conc_fresh[['RFd']] + man_pars_VS_cattle3.0$conc_fresh[['CP']] + 
man_pars_VS_cattle3.0$conc_fresh[['iNDF']]) * 0.42

for (i in c('starch', 'xa_dead','Cfat','RFd','CP','iNDF')){
  man_pars_VS_cattle3.0$conc_fresh[i] <- 0
}

# digestate
man_pars_VS_digestate3.0$conc_fresh['VSd'] <- (man_pars_VS_digestate3.0$conc_fresh[['starch']] + man_pars_VS_digestate3.0$conc_fresh[['xa_dead']] +
                                              man_pars_VS_digestate3.0$conc_fresh[['Cfat']] + man_pars_VS_digestate3.0$conc_fresh[['RFd']] + man_pars_VS_digestate3.0$conc_fresh[['CP']])
                                              
for (i in c('starch', 'xa_dead','Cfat','RFd','CP','iNDF')){
  man_pars_VS_digestate3.0$conc_fresh[i] <- 0
}


wthr_pars3.0 <- list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 10)
chem_pars3.0 <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561, RFd = 1/0.8444792, iNDF = 1/0.8444792, starch = 1/0.8444792, 
                                  Cfat = 1/0.3117844, CP = 1/0.6541602, VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, CO2_aer = 1/0.436, CO2_sr = 1/1.2, CO2_ureo = 1/1.57,
                                  N_CP = 1/0.1014, C_xa = 1/0.3753125, C_RFd = 1/0.376, C_iNDF = 1/0.358, N_xa = 1/0.08754375,
                                  C_starch = 1/0.377, C_Cfat = 1/0.265, C_CP = 1/0.359 , C_VFA = 1/0.374, C_VSd = 1/0.344, C_N_urea = 1/0.429,
                                  frac_CP_xa = 0.835))

arrh_pars_pig3.0 <- list(lnA = c(VSd_A = 31.3),
                         E_CH4 = c(VSd_A = 81000), 
                         A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                         E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                         R = 8.314,  
                         VS_CH4 = 6.67,
                         scale_alpha_opt = list(VSd = 1.19844, notVSd = 1.693),
                         kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

arrh_pars_cattle3.0 <- list(lnA = c(VSd_A = 31.2),
                         E_CH4 = c(VSd_A = 81000), 
                         A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 1.878889, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                         E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 9682.088, RFd = 81052, VSd = 81557, urea = 81559),  # CP is for GRA treatment in MILK project
                         R = 8.314,  
                         VS_CH4 = 6.67,
                         scale_alpha_opt = list(VSd = 0.815577, notVSd = 0.3737),
                         kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))


arrh_pars_digestate3.0 <- list(lnA = c(VSd_A = 31.2),
                            E_CH4 = c(VSd_A = 81000), 
                            A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 1.878889, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                            E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 9682.088, RFd = 81052, VSd = 81557, urea = 81559),  
                            R = 8.314,  
                            VS_CH4 = 6.67,
                            scale_alpha_opt = list(VSd = 0.554, notVSd = 0.6878),
                            kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

pig_pars3.0 <- c(grp_pars_pig3.0, arrh_pars_pig3.0, man_pars_pig3.0)
cattle_pars3.0 <- c(grp_pars_cattle3.0, arrh_pars_cattle3.0, man_pars_cattle3.0)
digestate_pars3.0 <- c(grp_pars_digestate3.0, arrh_pars_digestate3.0, man_pars_digestate3.0)

pig_parsVS3.0 <- c(grp_pars_VS_pig3.0, arrh_pars_pig3.0, man_pars_VS_pig3.0)
cattle_parsVS3.0 <- c(grp_pars_VS_cattle3.0, arrh_pars_cattle3.0, man_pars_VS_cattle3.0)
digestate_parsVS3.0 <- c(grp_pars_VS_digestate3.0, arrh_pars_digestate3.0, man_pars_VS_digestate3.0)


save(grp_pars_pig3.0, file = '../data/grp_pars_pig3.0.rda')
save(grp_pars_cattle3.0, file = '../data/grp_pars_cattle3.0.rda')
save(grp_pars_digestate3.0, file = '../data/grp_pars_digestate3.0.rda')
save(grp_pars_VS_pig3.0, file = '../data/grp_pars_VS_pig3.0.rda')
save(grp_pars_VS_cattle3.0, file = '../data/grp_pars_VS_cattle3.0.rda')
save(grp_pars_VS_digestate3.0, file = '../data/grp_pars_VS_digestate3.0.rda')
save(arrh_pars_pig3.0, file = '../data/arrh_pars_pig3.0.rda')
save(arrh_pars_cattle3.0, file = '../data/arrh_pars_cattle3.0.rda')
save(arrh_pars_digestate3.0, file = '../data/arrh_pars_digestate3.0.rda')

save(mic_pars3.0, file = '../data/mic_pars3.0.rda')
save(man_pars_pig3.0, file = '../data/man_pars_pig3.0.rda')
save(man_pars_VS_pig3.0, file = '../data/man_pars_VS_pig3.0.rda')
save(man_pars_cattle3.0, file = '../data/man_pars_cattle3.0.rda')
save(man_pars_VS_cattle3.0, file = '../data/man_pars_VS_cattle3.0.rda')
save(man_pars_digestate3.0, file = '../data/man_pars_digestate3.0.rda')
save(man_pars_VS_digestate3.0, file = '../data/man_pars_VS_digestate3.0.rda')

save(wthr_pars3.0, file = '../data/wthr_pars3.0.rda')
save(chem_pars3.0, file = '../data/chem_pars3.0.rda')
save(pig_pars3.0, file = '../data/pig_pars3.0.rda')
save(pig_parsVS3.0, file = '../data/pig_parsVS3.0.rda')
save(cattle_pars3.0, file = '../data/cattle_pars3.0.rda')
save(cattle_parsVS3.0, file = '../data/cattle_parsVS3.0.rda')
save(digestate_pars3.0, file = '../data/digestate_pars3.0.rda')
save(digestate_parsVS3.0, file = '../data/digestate_parsVS3.0.rda')


