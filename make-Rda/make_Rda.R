# Creates parameter objects

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

save(grp_pars2.0, file = '../data/grp_pars2.0.rda')
