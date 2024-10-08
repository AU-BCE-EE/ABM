coverfun <- function(cover, scale = 1){
  if(is.na(cover)) cover <- "none" 
  rd <- data.frame(none = 0, tent = 0.8, straw = 0.33, natural_crust = 0.45,
                   clay_pebbles = 0.41, floating_PVC = 0.16, biocover_porous_sheet = 0.66,
                   corrugated_sheets = 0.46, lid = 0.06, oil = 0.14, 
                   peat = 0.24, wood_chips = 0.53)
  EF_NH3_cover <- 1 - rd[[cover]]
  
  if(cover != 'none'){
    
    EF_NH3_cover <- EF_NH3_cover * scale
    
  }
  
  return(EF_NH3_cover)
}

# EFs from Sommer et al. 2022, "model for calcualting ammonia emission from stored animal liquid manure"