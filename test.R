

mng_pars = list(slurry_prod_rate = 5700,   # kg/d
                slurry_mass = 39000,           # Initial slurry mass (kg) NTS: convert to depth??
                storage_depth = 0.6,         # Storge structure depth, assued to be maximum slurry depth (m)
                resid_depth = 0.05,         # Residual slurry depth after emptying (m)
                floor_area = 650,           # Currently for NH3 loss from barn floor (nothing to do with pit/tank floor)
                area = 715,                 # Area (assume vertical sides, but slurry also underneath the walking path) (m2)
                empty_int = 42,            # (days, every 6th week)
                temp_C = 20,
                wash_water = 75000,            
                wash_int = NA,
                rest_d = 5,
                cover = NA,
                resid_enrich = 0.9,
                slopes = c(urea = NA, slurry_prod_rate = NA),
                graze = list(start = 'May', duration = 100, hours_day = 8),
                scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1))
out <- abm(365, mng_pars = mng_pars)#add_pars = list(graze.duration = 150, graze.hours_day = 12))
mng_pars$graze
