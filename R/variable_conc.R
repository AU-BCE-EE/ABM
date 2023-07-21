variable_conc <- function(conc_fresh, xa_fresh, t, t_run){

if (is.data.frame(conc_fresh)){
  names(conc_fresh) <- gsub(pattern = "conc_fresh.", replacement = "", x = names(conc_fresh))
  names <- names(conc_fresh[,!grepl("time", names(conc_fresh))])
  tconc_fresh <- t(conc_fresh)[-c(1),]
  ttime <- t(conc_fresh)["time",]
  ID <- NULL
  for (i in 1:(length(conc_fresh)-1)){
    conc_fresh_fun <- approxfun(ttime, tconc_fresh[i,], method = 'linear', f = 1,
                                yleft = tconc_fresh[i,1], yright = tconc_fresh[i,ncol(tconc_fresh)], rule = 2)
    ID1 <- conc_fresh_fun(t + t_run)
    ID <- rbind(ID,ID1)
  }
  conc_fresh <- as.list(ID)
  names(conc_fresh) <- names
}


if (is.data.frame(xa_fresh)){
  names(xa_fresh) <- gsub(pattern = "xa_fresh", replacement = "", x = names(xa_fresh))
  names <- names(xa_fresh[,!grepl("time", names(xa_fresh))])
  txa_fresh <- t(xa_fresh[, names])
  ttime <- t(xa_fresh)["time",]
  ID <- NULL
  for (i in 1:(length(xa_fresh)-1)){
    xa_fresh_fun <- approxfun(ttime, txa_fresh[i,], method = 'linear', f = 1,
                              yleft = txa_fresh[i,1], yright = txa_fresh[i, ncol(txa_fresh)], rule = 2)
    ID1 <- xa_fresh_fun(t + t_run)
    ID <- rbind(ID,ID1)
  }
  xa_fresh <- c(ID)
  names(xa_fresh) <- names
}
  
return(c(list(conc_fresh = conc_fresh, xa_fresh = xa_fresh)))
}

