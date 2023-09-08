makeConcFunc <- function(dat) {
  
  if(is.data.frame(dat)){
  
  funcs <- list()  # Initialize an empty list to store the functions
    
  names(dat) <- gsub(pattern = "conc_fresh.", replacement = "", x = names(dat))
  names <- names(dat[,!grepl("time", names(dat))])
  tconc_fresh <- t(dat)[-c(1), ]
  ttime <- t(dat)["time", ]
  
  for (i in 1:(length(dat) - 1)) {
    fun_name <- paste0('conc_fresh_fun_', names[i])
    fun <- approxfun(ttime, tconc_fresh[i, ], method = 'linear', f = 1,
                     yleft = tconc_fresh[i, 1], yright = tconc_fresh[i, ncol(tconc_fresh)], rule = 2)
    
    # Assign the function to the list using the function name as a list element name
    funcs[[fun_name]] <- fun
  }
  } else{
  
  funcs <- dat
  
  }
  
  return(funcs)  # Return the list of functions
}