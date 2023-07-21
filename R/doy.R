doy <- function(month, begin){
  doyd <- data.frame(month = c("January", "February", "March", "April",
                               "May", "June", "July", "August", "September",
                               "October", "November", "December"),
                     doy = c(1, 32, 60, 91, 121, 152, 183, 213, 244,
                             274, 305, 335))
 doys <- 0
 days <- 0
 
  for(i in 1:length(month)){
    doys[i] <- doyd$doy[doyd$month %in% month[i]] 
    if (doys[i] >= doyd$doy[doyd$month == begin]){
    days[i] <- doys[i] - doyd$doy[doyd$month == begin]
    } else {
    days[i] <- doys[i] + 365 - doyd$doy[doyd$month == begin]
    }
  }
 
  return(list(doy = doys, day = days))
}




