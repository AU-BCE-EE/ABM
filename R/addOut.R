# Add interval output to earlier results
# new can be lsoda output or state variable vector (after emptying)

addOut <- function(main, new = NULL) {
  
  # Change output from matrix to data frame
  # Do not drop first (time 0) row
  if (inherits(new, 'matrix')) {
    new <- data.frame(new)
  } 
  
  if (inherits(main, 'matrix')) {
    main <- data.frame(main)
  } 

  if (inherits(new, 'data.frame')) {
    # Get previous time from main, if it exists
    if ('time' %in% names(main)) {
      t_add <- max(main$time)
    } else {
      t_add <- 0
    }
    # Change time in output to cumulative time for complete simulation
    new$time <- new$time + t_add
    # Add results to earlier ones
    main <- rbind(main, new)
  } else if (inherits(new, 'numeric')){
    # Add row to main, duplicate last row
    main <- main[c(1:nrow(main), nrow(main)), ]
    # Replace values that are present in new
    main[nrow(main), names(new)] <-new 
  } else {
    stop('Class of new is not data frame or vector but is ', class(new))
  }
  
  return(main)

}
