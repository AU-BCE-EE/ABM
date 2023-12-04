\name{parameters}

\alias{parameters}
\alias{arrh_pars2.0}
\alias{chem_pars2.0}
\alias{grp_pars2.0}
\alias{man_pars2.0}
\alias{mic_pars2.0}
\alias{wthr_pars2.0}
\alias{arrh_pars1.0}
\alias{chem_pars1.0}
\alias{grp_pars1.0}
\alias{man_pars1.0}

\docType{data}
\title{
  Default parameter sets for ABM
}

\description{
  Default parameter sets for the \code{\link{abm}} function.
}

\usage{
  data(grp_pars2.0)
  data(grp_pars1.0)
  data(man_pars2.0)
  data(man_pars1.0)
  data(mic_pars2.0)
  data(wthr_pars2.0)
  data(chem_pars2.0)
  data(arrh_pars2.0)
  
}

\format{

  These are default parameter sets.
  Sets with `2.0` are the latest default sets.
  `1.0` versions are from ...
  
}

\source{
Dalby, F.R., Hafner, S.D., Petersen, S.O., Vanderzaag, A., Habtewold, J., Dunfield, K., Chantigny, M.H., Sommer, S.G., 2021. A mechanistic model of methane emission from animal slurry with a focus on microbial groups. PLOS ONE 16, e0252881. \doi{10.1371/journal.pone.0252881}

Dalby, F.R., Hansen, M.J., Guldberg, L.B., Hafner, S.D., Feilberg, A., 2023. Simple Management Changes Drastically Reduce Pig House Methane Emission in Combined Experimental and Modeling Study. Environ. Sci. Technol. \doi{10.1021/acs.est.2c08891}

Dalby, F.R., Ambrose, H.W.,  Poulsen, J.S., Nielsen, J.L., Adamsen, A.P.S., 2023. Pig slurry organic matter transformation and methanogenesis at ambient storage temperatures. JEQ. \doi{10.1002/jeq2.20512}
}

\examples{
  out1 <- abm(man_pars = man_pars2.0, grp_pars = grp_pars1.0)
}

\keyword{datasets}
