makeXaFreshFunc <- function(dat) {
  
  if(is.data.frame(dat)){
  
  funcs <- list()  # Initialize an empty list to store the functions
  
  names(dat) <- gsub(pattern = "xa_fresh", replacement = "", x = names(dat))
  names <- names(dat[,!grepl("time", names(dat))])
  txa_fresh <- t(dat[, names])
  ttime <- t(dat)["time",]
  
  for (i in 1:(length(dat) - 1)) {
    fun_name <- paste0('xa_fresh_fun_', names[i])
    fun <- approxfun(ttime, txa_fresh[i, ], method = 'linear', f = 1,
                     yleft = txa_fresh[i, 1], yright = txa_fresh[i, ncol(txa_fresh)], rule = 2)
    
    # Assign the function to the list using the function name as a list element name
    funcs[[fun_name]] <- fun
  }
  } else{
    
    funcs <- dat
    
  }
  
 
  return(funcs)  # Return the list of functions
}