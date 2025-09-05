rm(list =ls())

devtools::load_all('../../ABM')

slurry_mass_dat <- as.data.frame(fread('slurry_mass_dat.csv'))
temp_C_dat <- as.data.frame(fread('temp_C_dat.csv'))
pH_dat <- as.data.frame(fread('pH_dat.csv'))