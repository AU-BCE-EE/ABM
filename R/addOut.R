# Add interval output to earlier results

addOut <- function(dat, out, y.eff, grps, n_mic, t_run) {
  
    # Change output from matrix to data frame
    # Do not drop first (time 0) row
    out <- data.frame(out)

    # Add effluent results
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff

    # Change time in output to cumulative time for complete simulation
    out$time <- out$time + t_run
  
    # Add results to earlier ones
    dat <- rbind(dat, out)
    
    return(dat)

}
