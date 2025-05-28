
expandPars <- function(pars, elnms, parnms) {
  
  for (i in parnms) {
    ppo <- pars[[i]]
    p_nms <- names(pars[[i]])
    if (any(p_nms == 'default')) {
      pars[[i]][elnms] <- pars[[i]]['default']
      if (any(p_nms != 'default')) {
        pars[[i]][p_nms[p_nms != 'default']] <- ppo[p_nms[p_nms != 'default']]
      }
    }
    if (any(p_nms == 'all')) {
      pars[[i]][elnms] <- pars[[i]]['all']
    }
    # Fix order, drop default element if present, drop unused names
    pars[[i]] <- pars[[i]][elnms]
    # Check for missing values
    if (any(is.na(pars[[i]]))) stop('Missing grp_pars elements in ', i, '.')
  }

  return(pars)
  
}
  

