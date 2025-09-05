# Run benchmark scenarios with simpler1 branch

# Function
source('interpm.R')

# Packages
library(data.table)

# Variable data
pH_dat <- fread('speed_test_dat/pH_dat.csv')
slurry_mass_dat <- fread('speed_test_dat/slurry_mass_dat.csv')
temp_C_dat <- fread('speed_test_dat/temp_C_dat.csv')

# Round times as much as possible without getting duplicates
slurry_mass_dat[, time := round(time, 2)]

# Interpolate pH_dat and temp_C to mass times
var_dat <- slurry_mass_dat
var_dat[, pH := approx(pH_dat$time, pH_dat$pH, xout = time, rule = 2)$y] 
var_dat[, temp_C := approx(temp_C_dat$time, temp_C_dat$temp_C, xout = time, rule = 2)$y] 

# Input parameters
mng_pars <- list(storage_depth = 4,     
                 area = 100,              
                 resid_enrich = 1)

sub_pars <- list(subs = c('cellulose', 'protein', 'lipids'),
                 forms = c(cellulose = 'C6H10O5', protein = 'C4 H6.1 O1.2 N', 
                           lipids = 'C57 H104 O6', urea = 'CO(NH2)2'),
                 T_opt_hyd = c(default = 60),
                 T_min_hyd = c(default = 0),
                 T_max_hyd = c(default = 90),
                 hydrol_opt = c(lipids = 0.01, protein = 0.05, cellulose = 0.1),
                 sub_fresh = c(lipids = 3, protein = 20, cellulose = 35),
                 sub_init = c(lipids = 3, protein = 20, cellulose = 35))

grp_pars <- list(grps = c('m0', 'm1', 'm2'),
                 yield = c(default = 0.05),
                 xa_fresh = c(default = 0.05),
                 xa_init = c(default = 0.05),
                 dd_rate = c(default = 0.02),
                 ksv = c(default = 1),
                 qhat_opt = c(m0 = 1, m1 = 1, m2 = 2),
                 T_opt = c(m0 = 18, m1 = 18, m2 = 28),
                 T_min = c(m0 = 0, m1 = 6.41, m2 = 6.41),
                 T_max = c(m0 = 25, m1 = 25, m2 = 38))

man_pars <- list(VFA_fresh = c(CH3COOH = 2), pH = 7, dens = 1000)
chem_pars <- list(COD_conv = c(CH4 = 1/0.2507))

# Drop duplicate times

dim(var_dat)
dim(slurry_mass_dat)
var_pars <- list(var = var_dat)

# ABM
devtools::load_all()
system.time({
out1 <- abm(365, mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars, var_pars = var_pars)
})

png('mass1.png')
  plot(slurry_mass ~ time, type = 'l', data = out1)
dev.off()

png('pH1.png')
  plot(pH ~ time, type = 'l', data = out1)
dev.off()

png('CH41.png')
  plot(CH4_emis_rate ~ time, type = 'l', data = out1)
dev.off()

names(out1)

# Profiling for bottlenecks
Rprof('profile.out')
out1 <- abm(365, mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars, var_pars = var_pars)
Rprof(NULL)

sink('profile_summ.txt')
summaryRprof('profile.out')
sink()

# Add in speciation of TAN and H2S
man_pars2 <- list(comps = c('H2S', 'NH4p'),
                 comp_fresh = c(H2S = 0.01, NH4p = 2.5), 
                 VFA_fresh = c(CH3COOH = 2),
                 pH = 7, dens = 1000)

chem_pars2 <- list(COD_conv = c(CH4 = 1/0.2507),
                   specs = c('NH3', 'HSm', 'CH3COOm'),
                   mspec = c(NH3 = 'NH4p', HSm = 'H2S', CH3COOm = 'CH3COOH'),
                   lka = c(NH3 = '- 0.09046 - 2729.31/temp_K', 
                           HSm = '- 3448.7/temp_K + 47.479 - 7.5227 * log(temp_K)',
                           CH3COOm = '- 4.8288 + 21.42/temp_K'))

devtools::load_all()
system.time({
out2 <- abm(365, mng_pars = mng_pars, man_pars = man_pars2, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars2, var_pars = var_pars)
})

names(out2)

# Add NH3 emission
chem_pars3 <- list(COD_conv = c(CH4 = 1/0.2507),
                   specs = c('NH3', 'HSm', 'CH3COOm'),
                   mspec = c(NH3 = 'NH4p', HSm = 'H2S', CH3COOm = 'CH3COOH'),
                   lka = c(NH3 = '- 0.09046 - 2729.31/temp_K', 
                           HSm = '- 3448.7/temp_K + 47.479 - 7.5227 * log(temp_K)',
                           CH3COOm = '- 4.8288 + 21.42/temp_K'),
                   kl = c(NH3 = 0.0001) * 86400)

devtools::load_all()
system.time({
out3 <- abm(365, mng_pars = mng_pars, man_pars = man_pars2, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars3, var_pars = var_pars)
})

png('NH4p_conc.png')
plot(NH4p_conc ~ time, data = out2, type = 'l', ylim = c(0, 3))
lines(NH4p_conc ~ time, data = out3, col = 'red')
dev.off()

# Add CO2 emission
man_pars4 <- list(comps = c('H2S', 'NH4p', 'CO2'),
                 comp_fresh = c(H2S = 0.01, NH4p = 2.5, CO2 = 1), 
                 VFA_fresh = c(CH3COOH = 2),
                 pH = 7, dens = 1000)

chem_pars4 <- list(COD_conv = c(CH4 = 1/0.2507),
                   specs = c('NH3', 'HSm', 'CH3COOm', 'HCO3m'),
                   mspec = c(NH3 = 'NH4p', HSm = 'H2S', CH3COOm = 'CH3COOH', HCO3m = 'CO2'),
                   lka = c(NH3 = '- 0.09046 - 2729.31/temp_K', 
                           HSm = '- 3448.7/temp_K + 47.479 - 7.5227 * log(temp_K)',
                           CH3COOm = '- 4.8288 + 21.42/temp_K', 
			   HCO3m = '-2.778 + -353.5305 -0.06092*temp_K + 21834.37/temp_K + 126.8339*log10(temp_K) -1684915/temp_K^2'),
                   kl = c(NH3 = 0.0001, CO2 = 0.001) * 86400)

devtools::load_all()
system.time({
out4 <- abm(365, mng_pars = mng_pars, man_pars = man_pars4, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars4, var_pars = var_pars)
})

png('CO2_emis.png')
plot(CO2_emis_cum ~ time, data = out4, type = 'l')
dev.off()

png('CO2_CH4_ratio.png')
plot(CO2_emis_cum / CH4_emis_cum ~ time, data = out4, type = 'l')
dev.off()

# Profiling for bottlenecks
Rprof('profile4.out')
out4 <- abm(365, mng_pars = mng_pars, man_pars = man_pars4, grp_pars = grp_pars, sub_pars = sub_pars, chem_pars = chem_pars4, var_pars = var_pars)
Rprof(NULL)

sink('profile_summ4.txt')
summaryRprof('profile4.out')
sink()


