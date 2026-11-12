# Check that microbial group parameters have one element for each group, in the right order
# grp_par_nms = names of microbial group parameters

checkGrpNames <- function(pars, grp_par_nms) {

  for (i in grp_par_nms) {
    if (!identical(names(pars[[i]]), pars$grps)) {
      stop('Microbial groups, i.e., element names in `grp_pars`, must match `grps`. Problem with ', i, '.')
    }
  }

}
