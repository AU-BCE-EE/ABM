# Modified: 

molMass <- function(form, elements = NULL) {

  ## Check argument
  #checkArgClassValue(form, 'character')

  # Loop through all elements in form
  mmass <- NULL
  for(f in form) {
    # If and only if first letter of form is lowercase, entire string is capitalized
    if(grepl('^[a-z]', f)) f <- toupper(f) 

    # Get coefficients of formula
    fc <- readFormula(f)

    if (!is.null(elements)) {
      fc <- fc[intersect(names(fc), elements)]
    }

    # Check for unidentified element
    if(any(!names(fc) %in% names(atom.weights))) stop('One or more elements in \"form\" is not in the database. You can add it to the \"atom.weights\" vector if you want to modify the function code. Otherwise send a request to sasha.hafner@bce.au.dk.')

    # Calculate molar mass, using names of fc for indexing
    mmass <- c(mmass, sum(atom.weights[names(fc)]*fc))
  }

  #names(mmass) <- form

  return(mmass)
}
