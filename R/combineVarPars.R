
combineVarPars <- function(var_pars) {
  
  # If there is more than one data frame in var_pars, combine into one
  if (length(var_pars) > 1) {
    for (i in 2:length(var_pars)) {
      dv <- var_pars[[i]]
      nv <- names(var_pars)[[i]]
      if (nrow(dv) != nrow(var_pars$var)) {
        stop('var_pars has multiple elements with different sizes.')
      }
      if (!all(dv$time == var_pars$var$time)) {
        stop('Multiple elements in var_pars have different time values in first column.')
      }
      ll <- list()
      for (j in 1:nrow(var_pars$var)) {
        # It is challenging to get list into each element using indexing
        x <- as.numeric(dv[j, -1, drop = FALSE])
        names(x) <- names(dv)[-1]
        ll[[j]] <- x
      }
      var_pars$var[[nv]] <- ll
    }
  }

  # Take only the var element 
  # And var_pars must remain a list (see single brackets) to avoid duplicate par elements from combining into par
  var_pars <- var_pars['var']

  return(var_pars)

}
