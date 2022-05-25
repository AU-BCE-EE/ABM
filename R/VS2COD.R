# Sorts out fresh and initial concentrations based on TS/VS inputs

VS2COD <- function(x, cf) {

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

  # Convert to COD
  x[['COD']]   <- x[['VS']]   / cf
  x[['pCOD']]  <- x[['pVS']]  / cf
  x[['dpCOD']] <- x[['dpVS']] / cf
  x[['ipCOD']] <- x[['ipVS']] / cf
  x[['sCOD']]  <- x[['sVS']]  / cf
  x[['dsCOD']] <- x[['dsVS']] / cf
  x[['isCOD']] <- x[['isVS']] / cf

  return(x)

}
