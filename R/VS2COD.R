# Sorts out fresh and initial concentrations based on TS/VS inputs

VS2COD <- function(x, cf, sf = 0, digits = 6) {

  # cf = VS:COD conversion factor
  # sf = settling fraction

  # dsVS -> dsOM
  # dVSS -> dpOM

  # Rename a few 
  x[['dpVS']] <- x[['dVSS']] 
  x[['dsVS']] <- x[['dsVS']] 
  x[['pVS']] <- x[['VSS']] 

  # Calculate remainder by difference
  x[['sVS']]  <- x[['VS']] - x[['pVS']]
  x[['isVS']] <- x[['sVS']] - x[['dsVS']] 
  x[['ipVS']] <- x[['pVS']] - x[['dpVS']] 

  x[['iFS']] <- x[['TS']] - x[['VS']] 
  x[['ipFS']] <- x[['TSS']] - x[['VSS']] 
  x[['isFS']] <- x[['iFS']] - x[['ipFS']] 

  # Calculate settling component of degradable (sediment layer)
  x[['dpVSsed']] <- sf * x[['dpVS']] 
  x[['dpVS']] <- (1 - sf) * x[['dpVS']] 
  x[['ipVSsed']] <- sf * x[['ipVS']] 
  x[['ipVS']] <- (1 - sf) * x[['ipVS']] 
  x[['ipFSsed']] <- sf * x[['ipFS']] 
  x[['ipFS']] <- (1 - sf) * x[['ipFS']] 

  # Convert to COD
  x[['COD']]   <-    x[['VS']]      / cf
  x[['pCOD']]  <-    x[['pVS']]     / cf
  x[['dpCOD']] <-    x[['dpVS']]    / cf
  x[['dpCODsed']] <- x[['dpVSsed']] / cf
  x[['ipCOD']] <-    x[['ipVS']]    / cf
  x[['ipCODsed']] <- x[['ipVSsed']] / cf
  x[['sCOD']]  <-    x[['sVS']]     / cf
  x[['dsCOD']] <-    x[['dsVS']]    / cf
  x[['isCOD']] <-    x[['isVS']]    / cf

  x <- signif(x, digits)
  x[abs(x) < min(abs(x))/10^(digits - 3)] <- 0
  x[abs(x) < 1E-4] <- 0

  return(x)

}
