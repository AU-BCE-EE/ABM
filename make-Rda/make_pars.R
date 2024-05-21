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
                              ), 
                 kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

arrh_pars2.0 <- list(lnA = c(VSd_A = 31.3),
                     E_CH4 = c(VSd_A = 81000), 
                     A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                     E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                     R = 8.314,  
                     VS_CH4 = 6.67,
                     scale_alpha_opt = list(VSd = 1, notVSd = 1))

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
                    xa_fresh = c(default = 0.05628808),
                    xa_init = c(all = 0.05628808),
                    decay_rate = c(all = 0.02),
                    ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                    qhat_opt = c(m0 = 0.6903416, m1 = 0.2804, m2 = 2.014958, m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8.95),
                    T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                    T_min = c(m0 = 0, m1 = 13.59929, m2 = 13.59929, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                    T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                    ki_NH3_min = c(all = 0.015),
                    ki_NH3_max = c(all = 0.3),
                    ki_NH4_min = c(all = 2.7),
                    ki_NH4_max = c(all = 19.14),
                    ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                    ki_H2S_int = c(default = 1.03828, sr1 = 1.42/1.021 * 1.03828),
                    ki_H2S_min = c(default = 0.08),
                    IC50_low = c(default = 0.219629, sr1 = 0.2772/0.20854 * 0.219629),
                    ki_HAC = c(default = 0.22965),
                    pH_UL = c(default = 8.0),
                    pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_pig3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                       yield = c(default = 0.05, sr1 = 0.065),
                       xa_fresh = c(default = 0.05868828),
                       xa_init = c(all = 0.05868828),
                       decay_rate = c(all = 0.02),
                       ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                       qhat_opt = c(m0 = 0.7126616, m1 = 0.182, m2 = 1.978, m3 = 8.43, m4 = 8, m5 = 8, sr1 = 8.95),
                       T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                       T_min = c(m0 = 0, m1 = 12.726, m2 = 12.726, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                       T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                       ki_NH3_min = c(all = 0.015),
                       ki_NH3_max = c(all = 0.3),
                       ki_NH4_min = c(all = 2.7),
                       ki_NH4_max = c(all = 19.14),
                       ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                       ki_H2S_int = c(default = 1.0303, sr1 = 1.42/1.021 * 1.0303),
                       ki_H2S_min = c(default = 0.08),
                       IC50_low = c(default = 0.2478447, sr1 = 0.2772/0.20854 * 0.2478447),
                       ki_HAC = c(default = 0.147),
                       pH_UL = c(default = 8.0),
                       pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_VS_cattle3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                              yield = c(default = 0.05, sr1 = 0.065),
                              xa_fresh = c(default = 0.4529251),
                              xa_init = c(all = 0.4529251),
                              decay_rate = c(all = 0.02),
                              ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                              qhat_opt = c(m0 = 2.03072758, m1 = 0.0231, m2 = 1.1682497, m3 = 8.43, m4 = 8, m5 = 8, sr1 = 8.95),
                              T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                              T_min = c(m0 = 0, m1 = 13.88, m2 = 13.88, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                              T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                              ki_NH3_min = c(all = 0.015),
                              ki_NH3_max = c(all = 0.3),
                              ki_NH4_min = c(all = 2.7),
                              ki_NH4_max = c(all = 19.14),
                              ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                              ki_H2S_int = c(default = 1.021, sr1 = 1.42/1.021 * 1.021),
                              ki_H2S_min = c(default = 0.08),
                              IC50_low = c(default = 0.20854, sr1 = 0.2772),
                              ki_HAC = c(default = 0.31),
                              pH_UL = c(default = 8.0),
                              pH_LL = c(default = 6.5, sr1 = 5.5))

grp_pars_cattle3.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                        yield = c(default = 0.05, sr1 = 0.065),
                        xa_fresh = c(default = 0.13),
                        xa_init = c(all = 0.13),
                        decay_rate = c(all = 0.02),
                        ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                        qhat_opt = c(m0 = 3.818325, m1 = 0.01, m2 = 0.99214, m3 = 8.43, m4 = 8, m5 = 8, sr1 = 8.95),
                        T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                        T_min = c(m0 = 0, m1 = 5, m2 = 5, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                        T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                        ki_NH3_min = c(all = 0.015),
                        ki_NH3_max = c(all = 0.3),
                        ki_NH4_min = c(all = 2.7),
                        ki_NH4_max = c(all = 19.14),
                        ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                        ki_H2S_int = c(default = 0.9455, sr1 = 1.42/1.021 * 0.9455),
                        ki_H2S_min = c(default = 0.08),
                        IC50_low = c(default = 0.20854, sr1 = 0.2772),
                        ki_HAC = c(default = 0.31),
                        pH_UL = c(default = 8.0),
                        pH_LL = c(default = 6.5, sr1 = 5.5))



mic_pars3.0 <- list(ks_SO4 = 0.00694,
                    km_urea = 0.913,
                    decay_rate_xa = 0.02
                    )
                    

man_pars3.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.2, TAN = 0.0, starch = 5.25, 
                                      VFA = 1.7, xa_aer = 0, xa_bac = 0, xa_dead = 0, Cfat = 27.6, CP = 21.1, RFd = 25.4, iNDF = 11.3, VSd = 0, 
                                      VSd_A = 55, VSnd_A = 23.5, ash = 15), pH = 7, dens = 1000)

wthr_pars3.0 <- list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 10)

chem_pars3.0 <- list(COD_conv = c(CH4 = 1/0.2507, xa = 1/0.7069561, RFd = 1/0.8444792, iNDF = 1/0.8444792, starch = 1/0.8444792, 
                                  Cfat = 1/0.3117844, CP = 1/0.6541602, VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, CO2_aer = 1/0.436, CO2_sr = 1/1.2, CO2_ureo = 1/1.57,
                                  N_CP = 1/0.1014, C_xa = 1/0.3753125, C_RFd = 1/0.376, C_iNDF = 1/0.358, N_xa = 1/0.08754375,
                                  C_starch = 1/0.377, C_Cfat = 1/0.265, C_CP = 1/0.359 , C_VFA = 1/0.374, C_VSd = 1/0.344, C_N_urea = 1/0.429,
                                  frac_CP_xa = 0.835), 
kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

arrh_pars_pig3.0 <- list(lnA = c(VSd_A = 31.3),
                     E_CH4 = c(VSd_A = 81000), 
                     A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                     E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                     R = 8.314,  
                     VS_CH4 = 6.67,
                     scale_alpha_opt = list(VSd = 0.744457, notVSd = 1.563552183))

arrh_pars_pig3.0 <- list(lnA = c(VSd_A = 31.3),
                         E_CH4 = c(VSd_A = 81000), 
                         A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                         E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                         R = 8.314,  
                         VS_CH4 = 6.67,
                         scale_alpha_opt = list(VSd = 0.744457, notVSd = 1.563552183))

arrh_pars_cattle3.0 <- list(lnA = c(VSd_A = 31.3),
                         E_CH4 = c(VSd_A = 81000), 
                         A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                         E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                         R = 8.314,  
                         VS_CH4 = 6.67,
                         scale_alpha_opt = list(VSd = 0.185876, notVSd = 0.217138))

arrh_pars_cattle3.0 <- list(lnA = c(VSd_A = 31.3),
                            E_CH4 = c(VSd_A = 81000), 
                            A = c(xa_dead= 3.61383 * 10^12, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12, VSd = 3.61383 * 10^12, urea = 4.38*10^15), 
                            E = c(xa_dead= 81557, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052, VSd = 81557, urea = 81559),  
                            R = 8.314,  
                            VS_CH4 = 6.67,
                            scale_alpha_opt = list(VSd = 0.185876, notVSd = 0.217138))


pig_pars3.0 <- c(grp_pars_pig3.0, arrh_pars_pig3.0)
cattle_pars3.0 <- c(grp_pars_cattle3.0, arrh_pars_cattle3.0)
pig_parsVS3.0 <- c(grp_pars_VS_pig3.0, arrh_pars_pig3.0)
cattle_parsVS3.0 <- c(grp_pars_VS_cattle3.0, arrh_pars_cattle3.0)

save(grp_pars_pig3.0, file = '../data/grp_pars_pig3.0.rda')
save(grp_pars_cattle3.0, file = '../data/grp_pars_cattle3.0.rda')
save(grp_pars_VS_pig3.0, file = '../data/grp_pars_VS_pig3.0.rda')
save(grp_pars_VS_cattle3.0, file = '../data/grp_pars_VS_cattle3.0.rda')
save(arrh_pars_pig3.0, file = '../data/arrh_pars_pig3.0.rda')
save(arrh_pars_cattle3.0, file = '../data/arrh_pars_cattle3.0.rda')

save(mic_pars3.0, file = '../data/mic_pars3.0.rda')
save(man_pars3.0, file = '../data/man_pars3.0.rda')
save(wthr_pars3.0, file = '../data/wthr_pars3.0.rda')
save(chem_pars3.0, file = '../data/chem_pars3.0.rda')
save(pig_pars3.0, file = '../data/pig_pars3.0.rda')
save(pig_parsVS3.0, file = '../data/pig_parsVS3.0.rda')
save(cattle_pars3.0, file = '../data/cattle_pars3.0.rda')
save(cattle_parsVS3.0, file = '../data/cattle_parsVS3.0.rda')



