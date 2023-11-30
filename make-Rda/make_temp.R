library(dplyr)

# make data.frame for outside temperature

weather_dat_DK <- read.csv('monthly_summary.csv')
save(weather_dat_DK, file = '../data/weather_dat_DK.rda')

outside_slurry_temp_NIR <- weather_dat_DK %>% filter(decade >= 2000) %>% group_by(month) %>% 
  summarise(temp_air_C = mean(temp.mean)) %>% 
  mutate(temp_C = 0.511 * temp_air_C + 5.1886) %>% mutate(time = c(15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350))
outside_slurry_temp_NIR <- 
  as.data.frame(rbind(outside_slurry_temp_NIR, data.frame(month = c(0,13), 
                                            temp_air_C = rep(mean(outside_slurry_temp_NIR$temp_air_C[outside_slurry_temp_NIR$month %in% c(1, 12)]), 2),
                                            temp_C = rep(mean(outside_slurry_temp_NIR$temp_air_C[outside_slurry_temp_NIR$month %in% c(1, 12)]) * 0.511 + 5.1886, 2),
                                            time = c(0, 365))) %>% arrange(time)) %>% select(time, temp_C)

save(outside_slurry_temp_NIR, file = '../data/outside_slurry_temp_NIR.rda')

outside_slurry_temp_dig_NIR <- weather_dat_DK %>% filter(decade >= 2000) %>% group_by(month) %>% 
  summarise(temp_air_C = mean(temp.mean)) %>% 
  mutate(temp_C = 0.75 * temp_air_C + 6.23) %>% mutate(time = c(15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350))
outside_slurry_temp_dig_NIR <- 
 as.data.frame(rbind(outside_slurry_temp_dig_NIR, data.frame(month = c(0,13), 
                                            temp_air_C = rep(mean(outside_slurry_temp_dig_NIR$temp_air_C[outside_slurry_temp_dig_NIR$month %in% c(1, 12)]), 2),
                                            temp_C = rep(mean(outside_slurry_temp_dig_NIR$temp_air_C[outside_slurry_temp_dig_NIR$month %in% c(1, 12)]) * 0.75 + 6.23, 2),
                                            time = c(0, 365))) %>% arrange(time)) %>% select(time, temp_C)

save(outside_slurry_temp_dig_NIR, file = '../data/outside_slurry_temp_dig_NIR.rda')


outside_slurry_temp_vechi <- 
  data.frame(time = c(0, 15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350, 365),
             temp_C = c(7.7, 7.4, 7.2, 8.6, 11.9, 14.9, 17.3, 19.4, 19.2, 16.7, 13.4, 10.6, 8, 7.7))

save(outside_slurry_temp_vechi, file = '../data/outside_slurry_temp_vechi.rda')

outside_slurry_temp_dig_vechi <- 
  data.frame(time = c(0, 15, 45, 76, 106, 137, 167, 198, 228, 259, 289, 320, 350, 365),
             temp_C = c(18.75, 18, 16.5, 12.75, 9, 11, 18, 23, 27.5, 20, 18.67, 17, 19.5, 18.75))

save(outside_slurry_temp_dig_vechi, file = '../data/outside_slurry_temp_dig_vechi.rda')

