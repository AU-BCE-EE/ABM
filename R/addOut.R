# Add interval output to earlier results

addOut <- function(main, new, t_add, y.eff = NULL) {
  
    # Change output from matrix to data frame
    # Do not drop first (time 0) row
    new <- data.frame(new)

    # Add effluent results
    if (!is.null(y.eff)) {
      new[, names(y.eff)] <- 0
      new[nrow(new), names(y.eff)] <- y.eff
    }

    # Change time in output to cumulative time for complete simulation
    new$time <- new$time + t_add

    # Add results to earlier ones
    main <- rbind(main, new)
    
    return(main)

}
