library(dplyr)
# Creates parameter objects
# best fit grp_pars2.0 based on multiple datasets using CP, RFd, starch and Cfat as inputs
grp_pars2.0 <- list(grps = c('m0', 'm1', 'm2','sr1'),
                 yield = c(default = 0.05, sr1 = 0.065),
                 xa_fresh = c(default = 0.0628),
                 xa_init = c(all = 0.0628),
                 decay_rate = c(all = 0.02),
                 ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                 qhat_opt = c(m0 = 0.5515672, m1 = 0.5883362, m2 = 1.5692192, m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8),
                 T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                 T_min = c(m0 = 0, m1 = 8.1928624, m2 = 8.1928624, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                 T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                 ki_NH3_min = c(all = 0.015),
                 ki_NH3_max = c(all = 0.13),
                 ki_NH4_min = c(all = 2.7),
                 ki_NH4_max = c(all = 4.8),
                 ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                 ki_H2S_int = c(default = 0.93066, sr1 = 1.2938),
                 ki_H2S_min = c(default = 0.08))

mic_pars2.0 <- list(ks_SO4 = 0.00694,
                               km_urea = 0.913,
                               alpha_opt = c(urea = 60, VSd = 0.04954023),
                               alpha_T_min = c(urea = 0, VSd = 0),
                               alpha_T_opt = c(urea = 50, VSd = 50),
                               alpha_T_max = c(urea = 60, VSd = 60))

man_pars2.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.2, TAN = 0.0, starch = 5.25, 
                                  VFA = 1.7, xa_dead = 0, Cfat = 27.6, CP = 21.1, RFd = 25.4, iNDF = 11.3, VSd = 0, 
                                  VSd_A = 55, VSnd_A = 23.5, ash = 15), pH = 7, dens = 1000)

wthr_pars2.0 <- list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 10)

chem_pars2.0 <- list(COD_conv = c(CH4 = 1/0.2507, xa_dead = 1/0.73, RFd = 1/0.8444792, iNDF = 1/0.8444792, starch = 1/0.8444792, 
                              Cfat = 1/0.3117844, CP = 1/0.6541602, VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, CO2_aer = 1/0.436, CO2_sr = 1/1.2, CO2_ureo = 1/1.57,
                              N_CP = 1/0.1014, C_xa_dead = 1/0.358, C_RFd = 1/0.376, C_iNDF = 1/0.358,
                              C_starch = 1/0.377, C_Cfat = 1/0.265, C_CP = 1/0.359 , C_VFA = 1/0.374, C_VSd = 1/0.344, C_N_urea = 1/0.429), 
                 kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02))

arrh_pars2.0 <- list(lnA = c(VSd_A = 31.3),
                 E_CH4 = c(VSd_A = 81000), 
                 A = c(xa_dead= 8.56*10^7, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12), 
                 E = c(xa_dead= 60600, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052),  
                 R = 8.314,  
                 VS_CH4 = 6.67)

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
                               ki_H2S_int = c(default = 0.93066, sr1 = 1.2938),
                               ki_H2S_min = c(default = 0.08))

mic_pars1.0 <- list(ks_SO4 = 0.00694,
                    km_urea = 0.913,
                    alpha_opt = c(urea = 60, VSd = 0.04954023),
                    alpha_T_min = c(urea = 0, VSd = 0),
                    alpha_T_opt = c(urea = 50, VSd = 50),
                    alpha_T_max = c(urea = 60, VSd = 60))


man_pars1.0 <- list(conc_fresh = list(sulfide = 0.01, urea = 3.17, sulfate = 0.2, TAN = 0.0, starch = 0, 
                                      VFA = 2, xa_dead = 0, Cfat = 0, CP = 0, RFd = 0, iNDF = 0, VSd = 75, 
                                      VSd_A = 55, VSnd_A = 22, ash = 15), pH = 7, dens = 1000)

save(grp_pars1.0, file = '../data/grp_pars1.0.rda')
save(mic_pars1.0, file = '../data/mic_pars1.0.rda')
save(man_pars1.0, file = '../data/man_pars1.0.rda')

# make data.frame for outside temperature

weather_dat_DK <- read.csv('monthly_summary.csv')
save(weather_dat_DK, file = '../data/weather_dat_DK.rda')

outside_slurry_temp_NIR <- weather_dat_DK %>% filter(decade >= 2000) %>% group_by(month) %>% 
  summarise(air_temp = mean(temp.mean)) %>% 
  mutate(temp_C = 0.511 * air_temp + 5.1886) %>% mutate(time = c(15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350))
outside_slurry_temp_NIR <- 
  as.data.frame(rbind(outside_slurry_temp_NIR, data.frame(month = c(0,13), 
                                            air_temp = rep(mean(outside_slurry_temp_NIR$air_temp[outside_slurry_temp_NIR$month %in% c(1, 12)]), 2),
                                            temp_C = rep(mean(outside_slurry_temp_NIR$air_temp[outside_slurry_temp_NIR$month %in% c(1, 12)]) * 0.511 + 5.1886, 2),
                                            time = c(0, 365))) %>% arrange(time) %>% select(time, temp_C))

save(outside_slurry_temp_NIR, file = '../data/outside_slurry_temp_NIR.rda')

outside_slurry_temp_dig_NIR <- weather_dat_DK %>% filter(decade >= 2000) %>% group_by(month) %>% 
  summarise(air_temp = mean(temp.mean)) %>% 
  mutate(temp_C = 0.75 * air_temp + 6.23) %>% mutate(time = c(15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350))
outside_slurry_temp_dig_NIR <- 
 as.data.frame(rbind(outside_slurry_temp_dig_NIR, data.frame(month = c(0,13), 
                                            air_temp = rep(mean(outside_slurry_temp_dig_NIR$air_temp[outside_slurry_temp_dig_NIR$month %in% c(1, 12)]), 2),
                                            temp_C = rep(mean(outside_slurry_temp_dig_NIR$air_temp[outside_slurry_temp_dig_NIR$month %in% c(1, 12)]) * 0.75 + 6.23, 2),
                                            time = c(0, 365))) %>% arrange(time) %>% select(time, temp_C))

save(outside_slurry_temp_dig_NIR, file = '../data/outside_slurry_temp_dig_NIR.rda')


outside_slurry_temp_vechi <- 
  data.frame(time = c(0, 15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350, 365),
             temp_C = c(7.7, 7.4, 7.2, 8.6, 11.9, 14.9, 17.3, 19.4, 19.2, 16.7, 13.4, 10.6, 8, 7.7))

save(outside_slurry_temp_vechi, file = '../data/outside_slurry_temp_vechi.rda')

outside_slurry_temp_dig_vechi <- 
  data.frame(time = c(0, 15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350, 365),
             temp_C = c(18.75, 18, 16.5, 12.75, 9, 11, 18, 23, 27.5, 20, 18.67, 17, 19.5, 18.75))

save(outside_slurry_temp_dig_vechi, file = '../data/outside_slurry_temp_dig_vechi.rda')