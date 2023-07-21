graze_fun <- function(t, t_run, days, slurry_prod_rate, graze_int, graze_hours){
  
  y <- c(rep(graze_int, ceiling(days/365)))
  z <- rep(seq(from = 0, to = ceiling(days/365) - 1, by = 1), each = 2) * 365 + y
  mat <- matrix(z, nrow = ceiling(days/365), ncol = 2, byrow = T)
  
  if(any(mat[,1] < t + t_run & mat[,2] > t + t_run)){
    graze_factor <- 1 - graze_hours/24
    slurry_prod_rate <- graze_factor * slurry_prod_rate
  }
  
  return(slurry_prod_rate)

}

