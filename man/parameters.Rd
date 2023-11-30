\name{parameters}
\alias{parameters}
\docType{data}
\title{
  Default and alternative parameter sets for ABM
}

\description{
  Default parameter sets for \code{\link{abm}} function
}

\usage{
  data(grp_pars2.0), current default \code{grp_pars} default parameter set (v2.0) 
  data(grp_pars1.0), alternative \code{grp_pars} parameter set (v1.0) used in Dalby et al. (2023, see source section below)
  data(man_pars2.0), current default \code{man_pars} default parameter set which is based on an average Danish pig diet and with substrates defined as \code{RFd}, \code{CP}, \code{Cfat} and \code{starch} rather than \code{VSd}
  data(man_pars1.0), alternative \code{man_pars} parameter set (v1.0) based on \code{VSd} as substrate input
  data(mic_pars2.0), current default \code{mic_pars} parameter set (v2.0)
  data(wthr_pars2.0), current default \code{wthr_pars} parameter set (v2.0)
  data(chem_pars2.0), current default \code{chem_pars} parameter set (v2.0)
  data(arrh_pars2.0), current default \code{arrh_pars} parameter set (v2.0)
  
}

\format{

  
}

\source{
  Dalby et al. 2023. DOI: 10.1021/acs.est.2c08891
  
}

\examples{
  viewSpec(man_pars2.0)
  viewSpec(grp_pars2.0)
}

\keyword{datasets}
